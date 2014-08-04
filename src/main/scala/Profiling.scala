package hacone.AgIn

import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import java.io.File
import java.io.RandomAccessFile

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.ListBuffer
import scala.io.Source
import scala.math._
import scala.sys.process._

import Classifier._
import hacone.AgIn.IOManager._

object Profiling {

    // deprecated, only for compatibility check.
    // Make a IPD-profile around given position
    def pointProf(ipds: Array[PacBioIPD], idx: Int): Profile = {
      // if the position doesn't have enough neighborhood, generate zero-coverage profile
      val ZeroProf = new Profile( (-10 to 10).toList.map { i => (i, 0.0, 0, 0.0, 0)} )

      if (idx-10<=0 || idx+11>=ipds.length) return ZeroProf 

      val _prof = (-10 to 10).toList.map { i =>
        val (fi, ri) = (ipds(idx+i).fipd, ipds(idx+1-i).ripd)
        val (fc, rc) = (ipds(idx+i).fcov, ipds(idx+1-i).rcov)
        val (cover, mean) = (fc + rc, (fi * fc + ri * rc)/(fc + rc).toDouble)
        // handle NaN and (positive-)Infinity and winsorize
        if (mean.isNaN) (i, 0.0, 0, 0.0, 0) else (i, math.min(mean, 10.0), cover, 0.0, 2)
      }
      new Profile(_prof)
    }

    // Iterator version.
    // Make a IPD-profile around given position
    def makeProfile(ipds: ArrayIterator[PacBioIPD], idx: Int): Option[Profile] = {
      val _prof = for (i <- -10 to 10) yield {

        val (fi, fc) = if (0<=idx+i) {
          ipds(idx+i) match {
            case Some(ipd) => (ipd.fipd, ipd.fcov)
            // case None => { (0.0, 0); return None } // TODO: which is better ?
            case None =>  (0.0, 0) // TODO: for compatibility, this need to be like this
          }
        } else (0.0, 0)

        val (ri, rc) = if (0<=idx+1-i) {
          ipds(idx+1-i) match {
            case Some(ipd) => (ipd.ripd, ipd.rcov)
            case None =>  (0.0, 0)
          }
        } else (0.0, 0)

        val (cover, mean) = (fc + rc, (fi * fc + ri * rc)/(fc + rc).toDouble)
        if (mean.isNaN) (i, 0.0, 0, 0.0, 0) else (i, math.min(mean, 10.0), cover, 0.0, 2)
      }
      new Some(Profile(_prof.toList))
    }

  // when cannot decide, return false // TODO: not efficient
  def isCpG(sequence: String, pos: Int, isPlusStrand: Boolean): Boolean = {
    if (isPlusStrand) {
      if (pos <= 0 || sequence.length-1 <= pos) false
      else sequence(pos) == 'C' && sequence(pos+1) =='G'
    } else {
      if (pos <= 1 || sequence.length <= pos) false
      else sequence(pos-1) == 'C' && sequence(pos) == 'G'
    }
  }

  // find indices of all Cs of CpG in forward strand
  def findCpG(start: Int, end: Int, sequence: String): List[Int] = {
    (start until end).toList.filter { i => isCpG(sequence, i, true) }
  }

  def findAllCpG(sequence: String): List[Int] = findCpG(1, sequence.length, sequence)

  type LIP = (List[Int], List[Int])
  def partitionMethylationState(idxs: List[Int], bis: Array[Bisulfite]): LIP = {

    // println(bis.take(20).map(_.score).mkString(" "))

    // you can customize these

    def isMethyl(b: Bisulfite) = b.score >= 0.5
    def isUnMethyl(b: Bisulfite) = ! isMethyl(b)

    /*
    def isMethyl(b: Bisulfite) = b.score >= 0.8 
    def isUnMethyl(b: Bisulfite) = b.score <= 0.2 && b.score > 0
    */

    val methyl = ListBuffer.empty[Int]
    val unmethyl = ListBuffer.empty[Int]

    // TODO: idea of recursion is OK, let's clean up using case LT, GT ...
    @scala.annotation.tailrec
    def loop(is: List[Int], bs: Array[Bisulfite]): Unit = {
        if (is.isEmpty || bs.isEmpty) return
        else if (is.head == bs.head.position) {
          if (isMethyl(bs.head)) methyl += is.head
          else if (isUnMethyl(bs.head)) unmethyl += is.head
          loop(is.tail, bs.tail)
        }
        else if (is.head > bs.head.position) loop(is, bs.tail)
        else if (is.head < bs.head.position) loop(is.tail, bs)
    }

    loop(idxs, bis)

    (methyl.result, unmethyl.result)
  }

  // find indices of methylated Cs in forward strand
  def extractMCpG(idxs: List[Int], bis: Array[Bisulfite]): List[Int] = {
    def isMethylated(idx: Int): Boolean = {
      bis.find(_.position == idx) match {
        case Some(b) => if (b.score > 0.5) true else false
        case _ => false
      }
    }
    idxs.filter(isMethylated)
  }
  def extractUMCpG(idxs: List[Int], bis: Array[Bisulfite]): List[Int] = {
    def isUnMethylated(idx: Int): Boolean = {
      bis.find(_.position == idx) match {
        case Some(b) => if (b.score == 0) true else false
        case _ => false
      }
    }
    idxs.filter(isUnMethylated)
  }

  // take CpG position in forward strand, return all CpG in both strand
  def listBothStrandCpG(idxs: List[Int]): List[(Boolean, Int)] = {
    { for (i <- idxs) yield {
      List((true, i), (false, i+1))
    }}.foldLeft(List.empty[(Boolean, Int)])(_:::_)
  }

  // make profile around idxs, idxs: List[(isPlusStrand, index)] 
  // add filter : 10/28
  def profileOnIdxs(ipds: Array[PacBioIPD], idxs: List[(Boolean, Int)], profsize: Int, filt: (Double, Int) => Boolean = {(a,b) => true}): Profile = {
    val _prof = for (r <- -profsize to profsize) yield {
      var ipdsum = 0.0
      var covsum = 0
      var possum = 0
      for ((s, i) <- idxs) {
        if (s && i+r < ipds.length && i+r > 0) {
          if (filt(ipds(i+r).fipd, ipds(i+r).fcov)) { // filtering abnormal ipds
            ipdsum += ipds(i+r).fipd * ipds(i+r).fcov
            covsum += ipds(i+r).fcov
            possum += 1
          }
        } else if (i-r < ipds.length && i-r > 0) {
          if (filt(ipds(i-r).ripd, ipds(i-r).rcov)) {
            ipdsum += ipds(i-r).ripd * ipds(i-r).rcov
            covsum += ipds(i-r).rcov
            possum += 1
          }
        }
      }
      (r, ipdsum / covsum, covsum, 0.0, possum) // TODO : to imple var (use fvar, rvar)
    }
    new Profile(_prof.toList)
  }

  // make profile around all CpG on sequence
  def profileAll(ipds: Array[PacBioIPD], profsize: Int, sequence: String): Profile = {
     profileOnIdxs(ipds, listBothStrandCpG(findCpG(1, sequence.length, sequence)), profsize)
  }

  // make profile around all CpGs in specified Block
  def profileBlock(ipds: Array[PacBioIPD], profsize: Int, sequence: String, start: Int, end: Int, filt: (Double, Int)=>Boolean = (a,b)=>true): Profile = {
     profileOnIdxs(ipds, listBothStrandCpG(findCpG(start, end, sequence)), profsize, filt)
  }
  
  def profileUnMethyl(ipds: Array[PacBioIPD], bis: Array[Bisulfite], profsize: Int, sequence: String): Profile = {
   val a = findCpG(1, sequence.length, sequence) 
   val b = extractUMCpG(a, bis)
     println("#um %d".format(b.length))
   val c = listBothStrandCpG(b)
   profileOnIdxs(ipds, c, profsize)
  }

  def profileMethyl(ipds: Array[PacBioIPD], bis: Array[Bisulfite], profsize: Int, sequence: String): Profile = {
   val a = findCpG(1, sequence.length, sequence) 
   val b = extractMCpG(a, bis)
     println("#m %d".format(b.length))
   val c = listBothStrandCpG(b)
   profileOnIdxs(ipds, c, profsize)
  }

// TODO : very slow ?
  def weightedMeanAnswer(bis: Array[Bisulfite], start: Int, end: Int): (Double, Int) = {
    var scrsum = 0.0
    var covsum = 0
    for (b <- bis; if start <= b.position && b.position < end) {
      scrsum += b.score * b.coverage
      covsum += b.coverage
    }
    (scrsum, covsum)
  }
}
