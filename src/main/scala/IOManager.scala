package hacone.AgIn

import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import java.io.RandomAccessFile

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.ListBuffer
import scala.io.Source
import scala.math._

object IOManager {
  // read .input file (generated from csv)
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

  def readWigAsArray(fileScore: String, fileCover: String): Array[Bisulfite] = {
    val buffer = ArrayBuffer.empty[Bisulfite]
    for (((_score, _cover), _index) <- Source.fromFile(fileScore).getLines.zip(Source.fromFile(fileCover).getLines).zipWithIndex) {
      val score = _score.toDouble
      val cover = _cover.toInt
      val index = _index + 1 // zipWithIndex is 0-origined
      if (score != -1) buffer += Bisulfite(index, score, cover)
      // buffer += Bisulfite(index, score, cover)
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
}
