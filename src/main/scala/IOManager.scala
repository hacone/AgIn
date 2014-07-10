package hacone.AgIn

import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import java.io.RandomAccessFile

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.ListBuffer
import scala.io.Source
import scala.math._


object IOManager extends xerial.core.log.Logger {
  // read .input file (generated from csv) // deprecated !
  def readInputAsArray(filename: String): Array[PacBioIPD] = {
    val buffer = ArrayBuffer.empty[PacBioIPD]
    val reader = new CSVReader(new FileReader(filename))
    var last = (1, 0:Long, 0.0, 0) // strand, position, ipdratio, coverage
    val nullIPD = PacBioIPD(-1.0, 0, 0.0, -1.0, 0, 0.0)
    buffer += nullIPD // sentinel at position 0
    for (line <- Iterator.continually(reader.readNext).takeWhile(_ != null).map(_.toList)) {
      line match {
        case _str::_pos::_ipd::_cov::rest =>
          val str = _str.toInt
          val pos = _pos.toLong
          val ipd = _ipd.toDouble
          val cov = _cov.toInt

          if (pos > last._2+1) {
            if (last._1 == 0) buffer += PacBioIPD(last._3, last._4, 0.0, -1.0, 0, 0.0)
            for (_ <- last._2+1 until pos) buffer += nullIPD
            if (str == 1) buffer += PacBioIPD(-1.0, 0, 0.0, ipd, cov, 0.0)
          } else if (pos == last._2+1) {
            if (last._1 == 0) buffer += PacBioIPD(last._3, last._4, 0.0, -1.0, 0, 0.0)
            if (str==1) buffer += PacBioIPD(-1.0, 0, 0.0, ipd, cov, 0.0)
          } else {
            buffer += PacBioIPD(last._3, last._4, 0.0, ipd, cov, 0.0)
          }
          last = (str, pos, ipd, cov)
        case _ => ()
      }
    }
    buffer.toArray
  }

  // I want to use this like:  for ((refname, ipds) <- LoadIPDSets) blah blah ...
  type TaggedIPD = (String, Array[PacBioIPD])
  type IPDTuple = (Int, Long, Double, Int)
  // TODO: rewrite once more
  // TODO: should be rewritten using genuine iterator
  // def loadIPD(filename: String): Stream[TaggedIPD] = {
  import scala.collection.immutable.StreamIterator
  def loadIPD(filename: String): StreamIterator[TaggedIPD] = {

    def setupArray(linestr: Stream[List[String]]): Array[PacBioIPD] = {
      val nullIPD = PacBioIPD(-1.0, 0, 0.0, -1.0, 0, 0.0)

      // str are to be sorted as 10 10 10 ... 
      // but PacBioIPD are const'ed as (data for 0, data for 1)
      def takePair(linestr: Stream[List[String]]): (Long, PacBioIPD, Stream[List[String]])  = {
        def parseOne(line: List[String]): IPDTuple = line match {
          case _pos:: _str:: base:: score:: tMean:: tErr:: modelP:: _ipd:: _cov:: rest =>
            (_str.toInt, _pos.toLong, _ipd.toDouble, _cov.toInt)
          case _ => error("ill formed line parsing %s".format(filename)); (0,0,0.0,0)
        }
        val fst = parseOne(linestr.head)
        if (linestr.tail.isEmpty) {
          if (fst._1 == 0) (fst._2, PacBioIPD(fst._3, fst._4, 0.0, -1.0, 0, 0.0), linestr.tail)
          else (fst._2, PacBioIPD(-1.0, 0, 0.0, fst._3, fst._4, 0.0), linestr.tail)
        } else {
          val snd = parseOne(linestr.tail.head)
          if (fst._2 == snd._2) {
            if (fst._1 == 0) (fst._2, PacBioIPD(fst._3, fst._4, 0.0, snd._3, snd._4, 0.0),
                              linestr.tail.tail)
            else (fst._2, PacBioIPD(snd._3, snd._4, 0.0, fst._3, fst._4, 0.0),
                                    linestr.tail.tail)
          } else {
            if (fst._1 == 0) (fst._2, PacBioIPD(fst._3, fst._4, 0.0, -1.0, 0, 0.0), linestr.tail)
            else (fst._2, PacBioIPD(-1.0, 0, 0.0, fst._3, fst._4, 0.0), linestr.tail)
          }
        }
      }

      val buffer = ArrayBuffer.empty[PacBioIPD]
      @scala.annotation.tailrec
      def recur(n: Long, str: Stream[List[String]]): Unit = {
        if (str.isEmpty) return
        else {
          val (pos, ipd, next) = takePair(str)
          val diff: Int = (pos - n).toInt
          buffer ++= Array.fill(diff)(nullIPD)
          buffer += ipd
          recur(pos+1, next)
        }
      }

      recur(0, linestr)
      buffer.toArray
    }

    def body(linestr: Stream[List[String]]): Stream[TaggedIPD] = {
      if (linestr.isEmpty) Stream.empty[TaggedIPD]
      else {
        val _refname = linestr.head.head
        val refname = _refname.takeWhile(!_.isSpaceChar)

        info("setup IPD array for %s".format(refname))

        val (hd, tl) = linestr.span(_.head == _refname)

        info("spaned...")
        (refname, setupArray(hd.map(_.tail))) #:: body(tl)
      }
    }

    val reader = new CSVReader(new FileReader(filename))
    val ipdstr = body(
      Stream.continually(reader.readNext)
      .takeWhile(_ != null).map(_.toList).tail // chop header line
    )

    new scala.collection.immutable.StreamIterator(ipdstr)
  }

  def readWigAsArray(fileScore: String, fileCover: String): Array[Bisulfite] = {
    val buffer = ArrayBuffer.empty[Bisulfite]
    for (((_score, _cover), _index) <- Source.fromFile(fileScore).getLines.zip(
                                       Source.fromFile(fileCover).getLines)
                                       .zipWithIndex) {
      val (score, cover) = (_score.toDouble, _cover.toInt)
      val index = _index + 1 // adjust, since zipWithIndex is 0-origin
      if (score != -1.0) buffer += Bisulfite(index, score, cover)
    }
    buffer.toArray
  }

  def readWigWithConv(fileScore: String, fileCover: String, refname: String): Array[Bisulfite] = {
    val lot = {for (l <- Source.fromFile("input/LOT_nona_contig.dat").getLines.filter(_.contains(refname))) yield {
      l.split(' ').toList match {
        case n::os::oe::ns::ne::tp => (os.toInt, oe.toInt, ns.toInt, ne.toInt) // TODO: handle type
        case _ => (0, 0, 0, 0)
      }
    }}.toList.sortWith((x,y) => x._1 < y._1)

    val trim = {for (l <- Source.fromFile("input/LOT_trim.dat").getLines.filter(_.contains(refname))) yield {
      l.split(' ').toList match {
        case n::os::oe::tp => (os.toInt, oe.toInt) // TODO: handle type
        case _ => (0, 0)
      }
    }}.toList.sortWith((x,y) => x._1 < y._1)

    // println("lot len %d".format(lot.length))
    var wrong = 0

    def conv(oidx: Int): Int = {
      var overtrimmed = 0
      var ht = 0
      var tt = 0

      for (l <- lot) {
        if (oidx < l._1) return -1
        if (l._3 > l._4) {
          overtrimmed += l._3 - l._4
        }
        if (l._1 <= oidx && oidx < l._2) {
          for (t <- trim) {
            if (l._1 == t._1) ht = t._2 - t._1
            if (l._2 == t._2) tt = t._2 - t._1
          }
          if (l._3+(oidx-l._1-ht) != l._4-(l._2-oidx-tt)) wrong += 1
          return l._3 + (oidx - l._1 - ht) + overtrimmed
        }
      }
      return -1
    }

    val buffer = ArrayBuffer.empty[Bisulfite]
    var missing = 0
    for (((_score, _cover), _index) <- Source.fromFile(fileScore).getLines.zip(Source.fromFile(fileCover).getLines).zipWithIndex) {
      val score = _score.toDouble
      if (score != -1) {
        val cover = _cover.toInt
        val index = conv(_index) + 1 // zipWithIndex is 0-origined
        if (index != 0) {
          buffer += Bisulfite(index, score, cover)
        } else {
          missing += 1
        }
      }
    }
    // println("wrong = %d, missing = ".format(wrong, missing))
    buffer.toArray
  }

  // TODO : rewrite with Option
  // TODO : rename the member of tuple properly
  def readSequenceAsString(filename: String, refname: String): String = {
    // find index for refname from .fai file
    var refIndex = ("not found", 0, 0:Long, 0, 0)
    info("look for fasta: " ++ filename)
    Source.fromFile(filename + ".fai").getLines.find(_.contains(refname)) match {
      case Some(l) =>
        l.split('\t').toList match {
          // name, len, from, chars, bytes
          case n::l::f::c::b::rest => refIndex = (n, l.toInt, f.toLong, c.toInt, b.toInt)
          case _ => {
            refIndex = ("not found", 0, 0:Long, 0, 0)
            throw new Exception("bad format .fai")
          }
        }
      case _ => {
        refIndex = ("not found", 0, 0:Long, 0, 0)
        throw new Exception("not found .fai or column for specified refname")
      }
    }

    // read from fasta
    val buffer = new StringBuilder
    val rbuf = new Array[Byte](refIndex._5)
    val fasta = new RandomAccessFile(filename, "r")
    fasta.seek(refIndex._3)

    while (buffer.length < refIndex._2) {
      fasta.read(rbuf)
      buffer ++= new String(rbuf.take(refIndex._4))
    }

    "N" ++ buffer.toString.take(refIndex._2) // sentinel at index 0
  }

  import AgIn.Tsegment
  def writeGFF(outfile: String, refname: String, segments: List[Tsegment], command: String): Unit = {
    val gff = new java.io.PrintWriter(outfile)

    val AgInVersion = "0.9"

    // write out header
    gff.println("##gff-version 3")
    gff.println("##date %tc".format(new java.util.Date()))
    gff.println("##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12")
    gff.println("##source AgIn %s".format(AgInVersion)) // TODO: any other way to access this info ?
    gff.println("##source-commandline %s".format(command))

    for ((begin, end, avgscr, size) <- segments) {
      val line = List(
        refname, ".", "epigenetically_modified_region",
        begin.toString, end.toString,
        ".", ".", ".",
        "type=hypomethylated_interval;avg_coverage=%f;avg_score=%f;nCpG=%d;".format(-1.0, avgscr, size)
      ).mkString("\t")
      gff.println(line)
    }
  gff.close
  }
}
