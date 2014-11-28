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
        case _ => info("ill formed line: %s".format(line))
      }
    }
    // buffer.take(10000).foreach(x => println("%.2f %3d: %.2f %3d".format(x.fipd, x.fcov, x.ripd, x.rcov)))
    buffer.toArray
  }

  // I want to use this like:  for ((refname, ipds) <- LoadIPDSets) blah blah ...
  type TaggedIPD = (String, ArrayIterator[PacBioIPD])
  type IPDTuple = (Int, Long, Double, Int)

  case class ArrayIterator[T: reflect.ClassTag](val it: Iterator[T], val size: Int) {
    var hasNext = true
    var idx = 0 // next index of array to load into
    val arr = new Array[T](size)

    def apply(n: Int): Option[T] = {
      if (!(-1 < n && idx-size <= n)) {
        error("idx: %d size: %d n: %d".format(idx, size, n))
        assert(-1 < n && idx-size <= n)
      }
      while (n >= idx && it.hasNext) {
        arr(idx%size) = it.next()
        idx += 1
      }
      if (n < idx) {
        Some(arr(n%size))
      } else {
        hasNext = false
        None
      }
    }
  }

  // TODO (solved?): Exhaust remaining element in the Iterator
  // Split Iterator according to a predicate p (split when p gives false)
  class GroupIterator[T](it: Iterator[T], p: (T, T) => Boolean) extends Iterator[Iterator[T]] {
    val bit = it.buffered
    var last: Option[T] = None
    var first = true

    def hasNext(): Boolean =  bit.hasNext
    def next(): Iterator[T] = {
      // exhaust remaining element
      if (!first) {
      while (bit.hasNext && 
             { last match {
                 case Some(l) => p(l, bit.head)
                 case None => true
      }}) { last = Some(bit.next()) }
      } else { first = false }

      class Itr extends Iterator[T] {
        var _hasNext = true // TODO: you cant assume this !!
        def hasNext(): Boolean = _hasNext
        def next(): T = {
          val ret = bit.next()
          _hasNext = if (bit.hasNext) p(ret, bit.head) else false
          ret
        }
      }

      last = None
      if (bit.hasNext) new Itr else Iterator.empty
    }
  }

  // case class loadIPD(val filename: String) extends Iterator[TaggedIPD] {
  def loadIPD(filename: String): Iterator[TaggedIPD] = {

    type NameIPD = (String, IPDTuple)
    def parseOne(line: List[String]): NameIPD = line match {
      case rname::_pos:: _str:: base:: score:: tMean:: tErr:: modelP:: _ipd:: _cov:: rest =>
        (rname, (_str.toInt, _pos.toLong, _ipd.toDouble, _cov.toInt))
      case _ => {
        error("An ill-formed line when parsing \"%s\" in %s".format(line, filename));
        ("_ERROR_", (0,0,0.0,0))
      }
    }

    val reader = new CSVReader(new FileReader(filename))
    val it = {
      val _it = Iterator.continually(reader.readNext).takeWhile(_ != null)
      _it.next() // skip the header line that can't be parsed
      _it.map(x => parseOne(x.toList))
      .filter(_._1 != "_ERROR_").buffered
    }

    // TODO: this check should be treated in the constructor of GroupIterator
    for (_itForChr <- new GroupIterator(it, (a: NameIPD, b: NameIPD) => a._1 == b._1); itForChr = _itForChr.buffered; if itForChr.hasNext) yield {
      val ni = itForChr.head._2
      // info("next head: %d %d %.3f %d".format(ni._1, ni._2, ni._3, ni._4))
      (itForChr.head._1.takeWhile(!_.isWhitespace), setupArrIt(itForChr.map(_._2)))
    }

  }

    def setupArrIt(line: Iterator[IPDTuple]): ArrayIterator[PacBioIPD] = {
      class IPDIterator(it: Iterator[IPDTuple]) extends Iterator[PacBioIPD] {
        val nullIPD = PacBioIPD(-1.0, 0, 0.0, -1.0, 0, 0.0)
        var idx = 0 // next required position
        // var (buf1, buf2) = (nullIPD, nullIPD)
        var bufa: IPDTuple = (0,0,0.0,0)
        var bufb: IPDTuple = (0,0,0.0,0)

        var _hasNext = true
        def hasNext() = _hasNext

        def next(): PacBioIPD = {
          if (idx < bufa._2) {
            idx += 1
            nullIPD
          } else if (idx == bufa._2 && idx == bufb._2) {
            val ret = PacBioIPD(bufb._3, bufb._4, 0.0, bufa._3, bufa._4, 0.0)
            if (it.hasNext) bufa = it.next() else _hasNext = false
            if (it.hasNext) bufb = it.next() else bufb = (-1, -1, 0.0, 0)
            idx += 1
            ret
          } else if (idx == bufa._2) {
            val ret = if (bufa._1 == 0) {
              PacBioIPD(bufa._3, bufa._4, -1.0, 0, 0, 0.0)
            } else {
              PacBioIPD(-1.0, 0, 0.0, bufa._3, bufa._4, 0.0)
            }
            bufa = bufb
            if (bufa._1 == -1) _hasNext = false
            else if (it.hasNext) bufb = it.next()
            else bufb = (-1, -1, 0.0, 0)
            idx += 1
            ret
          } else throw new Exception("This could not happen")
        }
      }
      val result = new ArrayIterator(new IPDIterator(line), 100)
      result
    }

  // TODO: for now, wig must be in 1-origin
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

  // deprecated
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

    fasta.close

    "N" ++ buffer.toString.take(refIndex._2) // sentinel at index 0
  }

  import AgIn.Tsegment
  def writeHeaderToGFF(outfile: String, command: String): Unit = {
    val gff = new java.io.PrintWriter(new java.io.FileWriter(outfile+".gff", false))
    val AgInVersion = "0.9"

    // write out header
    gff.println("##gff-version 3")
    gff.println("##date %tc".format(new java.util.Date()))
    gff.println("##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12")
    gff.println("##source AgIn %s".format(AgInVersion)) // TODO: any other way to access this info ?
    gff.println("##source-commandline %s".format(command))
    gff.close
  }

  def writeRegionsToGFF(outfile: String, refname: String, segments: List[Tsegment]): Unit = {
    val gff = new java.io.PrintWriter(new java.io.FileWriter(outfile+".gff", true))

    for {
      (begin, end, avgscr, size) <- segments
      if avgscr < 0
    } {
      val line = List(
        refname, ".", "epigenetically_modified_region",
        begin.toString, end.toString,
        ".", ".", ".",
        // "type=hypomethylated_interval;avg_coverage=%f;avg_score=%f;nCpG=%d;".format(-1.0, avgscr, size)
        "type=hypomethylated_interval;avg_score=%.4f;nCpG=%d;".format(avgscr, size)
      ).mkString("\t")
      gff.println(line)
    }
    gff.close
  }

  // TODO: abstract them
  def writeCoverageToWig(outfile: String, refname: String,
                         ita: List[(Int, Double, Double)]): Unit = {
    val wig = new java.io.PrintWriter(new java.io.FileWriter(outfile+"_coverage.wig", true))
    wig.println("variableStep chrom=%s".format(refname))
    for ((idx, _, c) <- ita) {
      wig.println("%d %.2f".format(idx, c))
    }
    wig.close
  }

  def writeClassToWig(outfile: String, refname: String,
                         segment: List[(Int, Double, Int)]): Unit = {
    val wig = new java.io.PrintWriter(new java.io.FileWriter(outfile+"_class.wig", true))
    wig.println("variableStep chrom=%s".format(refname))
    for ((idx, _, c) <- segment) wig.println("%d %d".format(idx, c))
    wig.close
  }

  def writeContinuousPredictionToWig(outfile: String, refname: String,
                                     prediction: List[(Int, Int)]): Unit = {
    val wig = new java.io.PrintWriter(new java.io.FileWriter(outfile+"_continuous_class.wig", true))
    wig.println("variableStep chrom=%s".format(refname))
    for ((idx, c) <- prediction) wig.println("%d %d".format(idx, c))
    wig.close
  }
}
