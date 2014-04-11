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

  // when cannot decide, return false
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
    (for (i <- start until end; if isCpG(sequence, i, true)) yield i) toList
  }

  def findAllCpG(sequence: String): List[Int] = {
    (for (i <- 1 until sequence.length; if isCpG(sequence, i, true)) yield i) toList
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
