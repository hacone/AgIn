package hacone.AgIn

// TODO: this is messy, keep them close where they are used

import au.com.bytecode.opencsv.CSVReader
import java.io.FileReader
import java.io.File
import java.io.RandomAccessFile

import scala.collection.mutable.ArrayBuffer
import scala.collection.mutable.ListBuffer
import scala.collection.immutable.Map
import scala.io.Source
import scala.math._
import scala.sys.process._
import scala.util.parsing.json.JSON

import Classifier._
import hacone.AgIn.IOManager._
import hacone.AgIn.Profiling._
import hacone.AgIn.resources._

case class Bisulfite(position: Int, score: Double, coverage: Int)
case class PacBioIPD(fipd: Double, fcov: Int, fvar: Double, ripd: Double, rcov: Int, rvar: Double) // variation not yet used

object AgIn extends xerial.core.log.Logger {

  type Tsegment = (Int, Int, Double, Int) // begin, end, avg.score, size(#CpG)

  def main(args:Array[String]) {
    // here i place general settings used for every type of task
    // Some are set by convention, and others by configuration (config.dat file)

    /*
    val configures = Source.fromFile("./input/config.dat").getLines
                           .map(_.split('=').map(_.trim))
    val seglen = configures.find(_(0)=="L").getOrElse(Array("Default", "50"))(1).toInt
    val gamma = configures.find(_(0)=="Gamma").getOrElse(Array("Default", "-1.80"))(1).toDouble
    val fastapath = "./input/" + configures.find(_(0)=="Reference")
                                 .getOrElse(Array("Default", "GenRef.fasta"))(1)
    val vn = Resources.veclda
    */

    val ipdpath = "./input/IPD/"
    val wigcovpath = "./input/wigcov/"
    val wigscrpath = "./input/wigscr/"

    // calculate ita scores for each position from IPD profile, effective coverage, and Beta vector
    // TODO: short cut intermemiate result
    /*
    def prfstoIta(prfs: List[(Int, Profile)]): List[Double] = {
      def weightedProd(v: List[Double], u: List[Double], w: List[Double]): Double = {
        v.zip(u).zip(w).foldLeft(0.0)((s,t) => s + t._2*t._1._1*t._1._2)
      }
      {for(p <- prfs) yield weightedProd(
          {for (i <- -10 to 10) yield p._2(i) }.toList, vn,
          {for (i <- -10 to 10) yield p._2.cov(i).toDouble }.toList)
      }.toList
    }
    */

    // using ita, calculate final scores for given intercept. then call segmentation algorithm
    /*
    def callSegmentationIta(prfs: List[(Int, Profile)], ic: Double, ita: List[Double]): List[(Int, Double, Int)] = {
      val fm: Double => Double = { s => (s) }
      val fum: Double => Double = { s => -(s) } 
      val score = for ((i,p) <- ita.zip(prfs)) yield {
        i + 2 * p._2.mcv * ic
      }
      val optipdseg = OptInterval.OptInterval(score.toList, List(fum, fm), List(seglen,seglen))
      println("obtained optipdseg")

      return {for ((((a,x),b),c) <- prfs.zip(score).zip(optipdseg)) yield {
        (a, b, c) // position, score, class
      }}.toList
    }
    */

    // given segmentation, calculate sensitivity and other metrics, compared to bisulfite-seq as an answer set
    /*
    def evalPrediction(seg: List[(Int, Int)], bis: Array[Bisulfite]): Map[String, Double] = { // seg: List[(Position, class)]
      def sens(a: (Double,Double,Double,Double)): Double = { a._1/(a._1+a._4) }
      def prec(a: (Double,Double,Double,Double)): Double = { a._1/(a._1+a._3) }
      def spec(a: (Double,Double,Double,Double)): Double = { a._2/(a._2+a._3) }
      def mcc(a: (Double,Double,Double,Double)): Double = {
        (a._1 * a._2 - a._3 * a._4) /
        sqrt((a._1+a._3)*(a._1+a._4)*(a._2+a._3)*(a._2+a._4))
      }
      def f1(a: (Double,Double,Double,Double)): Double = { 2*prec(a)*sens(a)/(prec(a)+sens(a)) }
      var (tp,tn,fp,fn) = (0.0,0.0,0.0,0.0)
      var i = 0
      for (s <- seg) {
        while (bis(i).position < s._1) {
          i += 1
          if (i == bis.length) {
            val t = (tp,tn,fp,fn)
            return scala.collection.immutable.Map[String, Double](
              "tp" -> tp, "fp" -> fp, "tn" -> tn, "fn" -> fn, "sens" -> sens(t),
              "spec" -> spec(t), "prec" -> prec(t), "f1" -> f1(t), "mcc" -> mcc(t),
              "nipd" -> seg.length.toDouble, "neval" -> (tp+tn+fp+fn))
          }
        }
        if (bis(i).position == s._1) {
          if (bis(i).score >= 0) {
            // val pm = bis(i).score // probability where this position is methylated
            val pm = if (bis(i).score < 0.5) 0 else 1 // i rather use simple model
            if (s._2 == 1) { // M: negative
              tn += pm
              fn += 1-pm
            } else { // UM: positive
              tp += 1-pm
              fp += pm
            }
          } 
        }
      }
      val t = (tp,tn,fp,fn)
      scala.collection.immutable.Map[String, Double](
        "tp" -> tp, "fp" -> fp, "tn" -> tn, "fn" -> fn, "sens" -> sens(t),
        "spec" -> spec(t), "prec" -> prec(t), "f1" -> f1(t), "mcc" -> mcc(t),
        "nipd" -> seg.length.toDouble, "neval" -> (tp+tn+fp+fn))
    }
    */
    
    // given segmentation, zip each CpG site with bisulfite score on it.
    // When bisulfite score is missing, a negative number is assigned anyway
    /*
    def zipBis(seg: List[(Int, Int)], bis: Array[Bisulfite]): List[(Int, Int, Double)] = {
      val b = scala.collection.mutable.ListBuffer.empty[(Int, Int, Double)]
      var i = 0

      for (s <- seg) {
        if (i < bis.length) {
          while (i < bis.length && bis(i).position < s._1) { i += 1 }
          if (i == bis.length) { b += ((s._1, s._2, -3.0)) } // beyond bis
          else if (bis(i).position == s._1) { b += ((s._1, s._2, bis(i).score)) }
          else if (bis(i).position > s._1) { b += ((s._1, s._2, -2.0)) }
        } else {
          b += ((s._1, s._2, -3.0)) // beyond bis
        }
      }
      b.result
    }

    // (position, score, class)
    // TODO: rewrite
    def segmentSegments(seg: List[(Int, Double, Int)]): List[List[(Int, Double, Int)]] = {
      val b = scala.collection.mutable.ListBuffer.empty[List[(Int, Double, Int)]]
      var now = seg.head._3
      var bb = scala.collection.mutable.ListBuffer.empty[(Int, Double, Int)]
      bb += seg.head
      if (seg.tail == Nil) {
        b += bb.result
        return b.result
      }
      for (s <- seg.tail) {
        if (s._3 != now) {
          b += bb.result
          now = s._3
          bb = scala.collection.mutable.ListBuffer.empty[(Int, Double, Int)]
          bb += s
        } else {
          bb += s
        }
      }
      b += bb.result
      b.result
    }

    // get list of (start_idx, end_idx, average score) of each predicted region
    def scoreSegments(segs: List[List[(Int, Double, Int)]], fm: Double=>Double): List[(Int, Int, Double)] = {
      val b = scala.collection.mutable.ListBuffer.empty[(Int, Int, Double)]
      for (seg <- segs) {
        var score = 0.0
        for (s <- seg) { score += fm(s._2) }
        score /= seg.length
        b += ((seg.head._1, seg.last._1, score))
      }
      b.result
    }
    */

    // start of main: declaring handler of each task specified by task/**.json
    // TODO: it may be better to separate these task-specific functions. definitely...

    // write out IPDR profile and bisulfite score of each CpG site for estimation of a vector beta
    /*
    def extractProfile(protocol: Map[String, Any]): Unit = {
      (protocol("input"), protocol("output")) match {
        case (inputList: List[List[String]], outfile: String) => {
          for (input <- inputList) {
            input match {
              case List(refname, ipd, score, cover) => {
                println(refname)
                val sequence = readSequenceAsString(fastapath, refname)
                val bis = readWigAsArray(wigscrpath + score, wigcovpath + cover)
                val ipds = readInputAsArray(ipdpath + ipd)
                val cpgs = findAllCpG(sequence)

                // write prof vs bis-score for LDA in R
                val prfs: List[(Int, Profile)] = {for (i <- cpgs) yield { // (position of CpG, s(k))
                  (i, pointProf(ipds, i))
                }}.toList
                println("writing profiles on %s to %s".format(refname, outfile))

                // TODO: Anyway, don't rely this task. Use rather standard LDA vector supplied somewhere.
                hacone.AgIn.PrintData.printProfVsBis(prfs, bis, outfile)
              }
              case _ => println("invalid input file list: %s".format(input)) 
            }
          }
        }
        case _ => println("missing or invalid input(output) field"); sys.exit
      }
    }
    */

    /*
    def trainBeta(protocol: Map[String, Any]): Unit = {
      (protocol("input"), protocol("output")) match {
        case (inputList: List[List[String]], outfile: String) => {
          for (input <- inputList) {
            input match {
              case List(refname, ipd, score, cover) => {
                val logfile = new java.io.PrintWriter(outfile + refname + ".dat")
                println("loading %s".format(refname))
                val sequence = readSequenceAsString(fastapath, refname)
                val bis = readWigAsArray(wigscrpath + score, wigcovpath + cover)
                // Conversion was necessarry for using upgraded refs of PBJelly 
                // val bis = readWigWithConv(wigscrpath + score, wigcovpath + cover, refname)
                val ipds = readInputAsArray(ipdpath + ipd)
                val cpgs = findAllCpG(sequence)
                // prfs are (position of CpG, s(k))
                val prfs: List[(Int, Profile)] = cpgs.map { i => (i, pointProf(ipds, i)) }.filter(x => x._2.mcv > 0) // mcv is mean corverage of single strand
                val ita = prfstoIta(prfs)

                def icopt(iclist: List[Double], sig: Double): List[(Double, Double)] = {
                  return {for (ic <- iclist) yield {
                    val optseg = callSegmentationIta(prfs, ic, ita)
                    val toeval = {for((o,p) <- optseg.zip(prfs) if p._2.mcv > 0) yield (o._1, o._3)}.toList
                    val evp = evalPrediction(toeval, bis)
                    logfile.println("%.3f %.4f %.4f %.4f %.4f %.4f %d %.4f".format(
                      ic, evp("sens"), evp("spec"), evp("prec"), evp("f1"), evp("mcc"),
                      0, toeval.length.toDouble / cpgs.length))
                    (ic, evp("f1"))
                  }}.toList
                }
                import hacone.AgIn.GridSearcher._
                var icl = GridSearcher((-8, 5),20)
                val deltai = 0.01
                def mid(a: Double ,b: Double) = {(a+b) * 0.5}
                while (icl.max - icl.min > deltai) {
                  println("\n icl search : %.5f to %.5f".format(icl.min, icl.max))
                  icl = GridSearcher.next(icopt(icl, 0), 0.003)
                  println("optimizing (ic = %f)".format(mid(icl.max,icl.min)))
                }
                logfile.close
              }
              case _ => println("invalid input file list: %s".format(input)) 
            }
          }
        }
        case _ => println("missing or invalid input(output) field"); sys.exit
      }
    }
    */

    /*
    def pointwisePrediction(protocol: Map[String, Any]): Unit = {
      (protocol("input"), protocol("output")) match {
        case (inputList: List[List[String]], outfile: String) => {
          for (input <- inputList) {
            input match {
              case List(refname, ipd, score, cover) => {
                val logfile = new java.io.PrintWriter(outfile + refname + ".dat")
                println("loading %s".format(refname))
                val sequence = readSequenceAsString(fastapath, refname)
                val bis = readWigAsArray(wigscrpath + score, wigcovpath + cover)
                // val bis = readWigWithConv(wigscrpath + score, wigcovpath + cover, refname)
                val ipds = readInputAsArray(ipdpath + ipd)
                val cpgs = findAllCpG(sequence)
                val prfs: List[(Int, Profile)] = cpgs.map { i => (i, pointProf(ipds, i)) }.filter(x => x._2.mcv > 0) // mcv is mean corverage of single strand
                val ita = prfstoIta(prfs)

                def icopt(iclist: List[Double], sig: Double): List[(Double, Double)] = {
                  return {for (ic <- iclist) yield {
                    val optseg = for ((i,p) <- ita.zip(prfs)) yield {
                      (p._1, i + 2 * p._2.mcv * ic,
                        if(i + 2 * p._2.mcv * ic > 0) 1 else 0)
                    }
                    //val optseg = callSegmentationIta(prfs, ic, ita)
                    val evp = evalPrediction({for (((c,d,e)) <- optseg) yield (c, e)}.toList, bis)
                    logfile.println("%.3f %.4f %.4f %.4f %.4f %.4f".format(ic, evp("sens"), evp("spec"), evp("prec"), evp("f1"), evp("mcc")))
                    (ic, evp("f1"))
                  }}.toList
                }
                import hacone.AgIn.GridSearcher._
                var icl = GridSearcher((-8, 5),20)
                val deltai = 0.01
                def mid(a: Double ,b: Double) = {(a+b) * 0.5}
                while (icl.max - icl.min > deltai) {
                  println("\n icl search : %.5f to %.5f".format(icl.min, icl.max))
                  icl = GridSearcher.next(icopt(icl, 0), 0.003)
                  println("optimizing (ic = %f)".format(mid(icl.max,icl.min)))
                }
                logfile.close
              }
              case _ => println("invalid input file list: %s".format(input)) 
            }
          }
        }
        case _ => println("missing or invalid input(output) field"); sys.exit
      }
    }
    */
    
    /*
    def getGCrate(sequence: String): List[(Int, Double)] = {
      var rate = 0.0
      var nrate = 0.0
      val a = new Array[Int](100)
      val n = new Array[Int](100)
      val b = new scala.collection.mutable.ListBuffer[(Int, Double)]
      var i = 0
      for (s <- sequence) {
        a(i%100) = if (s=='G' || s=='C') 1 else 0 
        n(i%100) = if (s=='N') 1 else 0 
        rate += ( a(i%100) - a((i+1)%100) )
        nrate += ( n(i%100) - n((i+1)%100) )
        if (i%100 == 0){
          if (nrate > 50) {
            b += ((i-50, -50))
          } else {
            b += ((i-50, rate/(100-nrate)))
          }
        }
        i += 1
      }
      b.result
    }

    def binCpGs(cpgs: List[Int], gc: List[(Int, Double)]): List[(Int, List[Int])] = {
      val a = new scala.collection.mutable.ListBuffer[(Int, ListBuffer[Int])]
      for (thre <- 5 to 100 by 5) a += ((thre, new scala.collection.mutable.ListBuffer[Int]))
      var rest = gc
      for (s <- cpgs) {
        while (rest.head._1 < s - 50 && rest.tail != Nil) rest = rest.tail
        var target = a
        if (rest.head._2 >= 0) { 
          for (_ <- 1 to (rest.head._2 * 20).toInt) target = target.tail
          target.head._2 += s
        }
      }
      a.result.map(x => { (x._1, x._2.result) })
    }

    def binAny[T](cpgs: List[(T, Int)], gc: List[(Int, Double)]): List[(Int, List[(T, Int)])] = {
      val a = new scala.collection.mutable.ListBuffer[(Int, ListBuffer[(T,Int)])]
      for (thre <- 5 to 100 by 5) a += ((thre, new scala.collection.mutable.ListBuffer[(T,Int)]))
      var rest = gc
      for (s <- cpgs) {
        while (rest.head._1 < s._2 - 50 && rest.tail != Nil) rest = rest.tail
        if (rest.head._2 >= 0) {
          a(min((rest.head._2*20).toInt, 19))._2 += s
        }
      }
      a.result.map(x => { (x._1, x._2.result) })
    }

    // Deprecated : These 2 function were defined to filter out new sequences added by PBJelly
    // for fair comparison of CpG coverage of bisulfite-seq and PacBio
    def filterNewseq[T](anys: List[(T, Int)], refname: String): List[(T, Int)] = { // filter out newly filled seq
      val lot = {for (l <- Source.fromFile("input/LOT_newseq.dat").getLines.filter(_.contains(refname))) yield {
      l.split(' ').toList match {
        case n::os::oe::ns::ne::tp => (ns.toInt, ne.toInt)
        case _ => (-1, -1)
      }}}.toList.filter(x=>x._2>0)sortWith((x,y) => x._1 < y._1)
      anys.filter(x => ! lot.exists(s => s._1 <= x._2 && x._2 < s._2))
    }

    def filterNewseqRev[T](anys: List[(T, Int)], refname: String): List[(T, Int)] = { // return newly filled seq
      val lot = {for (l <- Source.fromFile("input/LOT_newseq.dat").getLines.filter(_.contains(refname))) yield {
      l.split(' ').toList match {
        case n::os::oe::ns::ne::tp => (ns.toInt, ne.toInt)
        case _ => (-1, -1)
      }}}.toList.filter(x=>x._2>0)sortWith((x,y) => x._1 < y._1)
      anys.filter(x => lot.exists(s => s._1 <= x._2 && x._2 < s._2)) // inverse !!!
    }
    */

    /*
    def reportCoverage(protocol: Map[String, Any]): Unit = { // evaluation of coverage (comparison to bisulfite data)
      (protocol("input"), protocol("output")) match {
        case (inputList: List[List[String]], outfile: String) => {
          val logfile = new java.io.PrintWriter(outfile)

          val covbinpac = (1 to 20).toList.map(x => new Array[Int](51))
          val covbinbis = (1 to 20).toList.map(x => new Array[Int](51))

          for (input <- inputList) {
            input match {
              case List(refname, ipd, score, cover) => {
                val sequence = readSequenceAsString(fastapath, refname)
                val bis = readWigAsArray(wigscrpath + score, wigcovpath + cover).filter(_.coverage > 0)
                val ipds = readInputAsArray(ipdpath + ipd)
                val cpgs = findAllCpG(sequence)
                val gc = getGCrate(sequence)
                val prfs: List[(Int, Profile)] = cpgs.map(i => (i, pointProf(ipds, i))).filter(_._2.mcv > 0)

                // add 13/04/24
                // val fcpgs = filterNewseq(cpgs.map(x=>(x,x)), refname).map(x=>x._1)
                // println("%s %d".format(refname, fcpgs.length))

                val binnedCpG = binCpGs(cpgs, gc)
                val binnedPac = binAny(prfs.map(x=>(x._2,x._1)), gc)

                print(refname)

                for ( ((bin, pac), (cbin, cpg)) <- binnedPac.zip(binnedCpG)) {
                  val covsum = pac.foldLeft(0.0)((x,y)=>x + y._1.mcv * 2)
                  print(" %3d %6d %5d %5d %2.2f /".format(cbin, covsum.toInt, pac.length, cpg.length, covsum / cpg.length))
                  for (p <- pac) {
                    covbinpac((bin/5)-1)(min((p._1.mcv*2).toInt, 50)) += 1
                  }
                  covbinpac((bin/5)-1)(0) += cpg.length - pac.length // # of CpG sites on which IPD info was missing
                }

                val binnedBis = binAny(bis.map(b=>(b, b.position)).toList , gc)
                for ( ((bin, bis), (cbin, cpg)) <- binnedBis.zip(binnedCpG) ) {
                  val covsum = bis.foldLeft(0.0)((x,y)=>x + y._1.coverage)
                  print(" %3d %6d %5d %5d %2.2f /".format(cbin, covsum.toInt, bis.length, cpg.length, covsum / cpg.length))
                  for (b <- bis) {
                    covbinbis((bin/5)-1)(min(b._1.coverage, 50)) += 1
                  }
                  covbinbis((bin/5)-1)(0) += cpg.length - bis.length // # of CpG sites on which bisulfite reads were missing
                }
              }
              case _ => println("invalid input file list: %s".format(input)) 
            }
          }

          // print each bin vs cover
          for (i <- 0 to 50) { print("%2d ".format(i)) }
          println("")
          for (b <- 0 to 19) {
            print("%d".format((b+1)*5))
            for (i <- 0 to 50) { print("%7d ".format(covbinpac(b)(i))) }
            println("")
          }
          for (b <- 0 to 19) {
            print("%d".format((b+1)*5))
            for (i <- 0 to 50) { print("%7d ".format(covbinbis(b)(i))) }
            println("")
          }
          logfile.close
        }
        case _ => println("missing or invalid input(output) field"); sys.exit
      }
    }
    */

    /*
    def testOnUnseen(protocol: Map[String, Any]): Unit = { // evaluation of predictive power within unseen data (comparison to bisulfite data)
      (protocol("input"), protocol("output")) match {
        case (inputList: List[List[String]], outfile: String) => {
          val logfile = new java.io.PrintWriter(outfile)

          val covbinpac = (1 to 20).toList.map(x => new Array[Int](51))
          val covbinbis = (1 to 20).toList.map(x => new Array[Int](51))

          // logfile.println("scaf len cpglen paclen bislen evallen tp tn fp fn sens spec prec f1") // schema

          for (input <- inputList) {
            input match {
              case List(refname, ipd, score, cover) => {
                // println("loading %s".format(refname))
                val sequence = readSequenceAsString(fastapath, refname)
                // deprecated 
                // val bis = readWigWithConv(wigscrpath + score, wigcovpath + cover, refname).sortWith((a,b) => a.position < b.position)
                val bis = readWigAsArray(wigscrpath + score, wigcovpath + cover)
                val ipds = readInputAsArray(ipdpath + ipd)
                val cpgs = findAllCpG(sequence)
                val gc = getGCrate(sequence)
                val binnedCpG = binCpGs(cpgs, gc)

                // val prfs: List[(Int, Profile)] = cpgs.map(i => (i, pointProf(ipds, i))).filter(_._2.mcv > 2.5)
                val prfs: List[(Int, Profile)] = cpgs.map(i => (i, pointProf(ipds, i))).filter(_._2.mcv > 0)

                val binnedPac = binAny(prfs.map(x=>(x._2,x._1)), gc)
                print(refname)
                for ( ((bin, pac), (cbin, cpg)) <- binnedPac.zip(binnedCpG) ) {
                  val covsum = pac.foldLeft(0.0)((x,y)=>x + y._1.mcv * 2)
                  print(" %3d %6d %5d %5d %2.2f /".format(cbin, covsum.toInt, pac.length, cpg.length, covsum / cpg.length))
                  for (p <- pac) {
                    covbinpac((bin/5)-1)(min((p._1.mcv*2).toInt, 50)) += 1
                  }
                  covbinpac((bin/5)-1)(0) += cpg.length - pac.length // # of CpG sites on which IPD info was missing
                }

                // val binnedBis = binAny(bis.map(b=>(b, b.position)).toList , gc)
                // val binnedBis = binAny(bis.map(b=>(b, b.position)).toList.filter(_._1.coverage > 2.5) , gc)
                val binnedBis = binAny(bis.map(b=>(b, b.position)).toList.filter(_._1.coverage > 0) , gc)
                for ( ((bin, bis), (cbin, cpg)) <- binnedBis.zip(binnedCpG) ) {
                  val covsum = bis.foldLeft(0.0)((x,y)=>x + y._1.coverage)
                  print(" %3d %6d %5d %5d %2.2f /".format(cbin, covsum.toInt, bis.length, cpg.length, covsum / cpg.length))
                  for (b <- bis) {
                    covbinbis((bin/5)-1)(min(b._1.coverage, 50)) += 1
                  }
                  covbinbis((bin/5)-1)(0) += cpg.length - bis.length // # of CpG sites on which bisulfite reads were missing
                }

                // val ita = prfstoIta(prfs)

                // val optseg = callSegmentationIta(prfs, gamma, ita)
                // val scores = scoreSegments(segmentSegments(optseg), ((a:Double) => a))

                // TODO: TODO TODO TODO TODO TODO 
                for (ksi <- (0*20 to 10*20).map(_/40.0) ++
                               (10*10+1 to 20*10).map(_/20.0)) {
                  val toeval0 = {for (l <- scores.zip(segmentSegments(optseg)); if abs(l._1._3) > ksi) yield {
                    l._2
                  }}.flatten
                  val toeval = {for(o <- toeval0) yield (o._1, o._3)}.toList

                  val evp = evalPrediction(toeval, bis)
                  println("%s %d %d %d %d %d %d %d %d %d %.4f %.4f %.4f %.4f %.4f".format(
                    refname, sequence.length,
                    cpgs.length, prfs.length, bis.length, evp("neval").toInt,
                    evp("tp").toInt, evp("tn").toInt, evp("fp").toInt, evp("fn").toInt,
                    evp("sens"), evp("spec"), evp("prec"), evp("f1"), minabs
                  ))
                }

                for (nd <- 0 to 10) {
                val toeval = {
                  val ksi = 0.0 // should not be here ??
                  // val nd = 5
                  val goodSegments =
                  scores.zip(segmentSegments(optseg)).filter(x=>abs(x._1._3)>ksi).map(x=>x._2.drop(nd).dropRight(nd)).flatten
                  goodSegments.map(x => (x._1, x._3))
                }

                val evp = evalPrediction(toeval, bis)
                println("%d %s %d %d %d %d %d %d %d %d %d %.4f %.4f %.4f %.4f".format(nd,
                  refname, sequence.length,
                  cpgs.length, prfs.length, bis.length, evp("neval").toInt,
                  evp("tp").toInt, evp("tn").toInt, evp("fp").toInt, evp("fn").toInt,
                    evp("sens"), evp("spec"), evp("prec"), evp("f1")))
                }
                println("")

                // val filledhypo = findFilledHypo(segmentSegments(optseg), bis)
                // val zb = zipBis(toeval, bis)

                // write out result
                
                for ((((a,b,c),p),z) <- optseg.zip(prfs).zip(zb)) logfile.println("%d %f %d %f %f".format(a, b, c, p._2.mcv, z._3))
                val confile = new java.io.PrintWriter("IPDConfval.dat")
                for ((a, b, c, h) <- filledhypo) confile.println("%d %d %f %f".format(a, b, c, h))
                confile.close

                val f1: Double => Double = {s: Double => (s)}
                val f2: Double => Double = {s: Double => (-s)}

                for (i <- List(0.7)) yield {
                  val theta1 = i
                  val theta2 = 1-i

                  val f1: Bisulfite => Double = { b: Bisulfite =>
                    b.coverage * (b.score  - (1-b.score))
                  }
                  val f2: Bisulfite => Double = { b: Bisulfite =>
                    -b.coverage * (b.score - (1-b.score))
                  }

                  val d: List[(Int, Bisulfite)] = {for (b <- bis) yield {
                    if (b.score < 0) (-1, new Bisulfite(-1, -1.0, 0))
                    else (b.position, b)
                }}.toList.filter(_._1 > 0)

                  val optseg = OptInterval.OptInterval({for (dd<-d) yield dd._2}.toList, List(f2, f1), List(seglen, seglen))
                  assert(d.length == optseg.length, "opt seg failed ?")
                  
                  val evp = evalPrediction({for ((c,(d,e)) <- optseg.zip(d)) yield (d, c)}.toList, bis)
                  println("answer performance\n# %.4f %.4f %.4f %.4f %.4f".format(evp("sens"), evp("spec"), evp("prec"), evp("f1"), evp("mcc")))

                  val anssegfile = new java.io.PrintWriter("BisulfiteSegment%d.dat".format((i*10).toInt))
                  val ansseg = {for ((a, b) <- d.zip(optseg)) yield {
                    anssegfile.println("%d %f %d".format(a._1, a._2.score, b))
                    (a._1, a._2.score, b)
                  }}.toList
                  anssegfile.close

                  val gcfile = new java.io.PrintWriter("GCrate.dat")
                  for (r <- gc) { gcfile.println("%d %f".format(r._1, r._2)) }
                  gcfile.close
                }
              }
              case _ => println("invalid input file list: %s".format(input)) 
            }
          }
          logfile.close
        }
        case _ => println("missing or invalid input(output) field"); sys.exit
      }
    }
    */

    /*
    def predictOnRealCase(protocol: Map[String, Any]): Unit = {
      info("gamma = %f".format(gamma))
      (protocol("input"), protocol("output")) match {
        case (inputList: List[List[String]], outfile: String) => {
          inputList.map { input: List[String] =>
            input match {
              case List(refname, ipd) => {
                val logfile = new java.io.PrintWriter(outfile + refname + ".allpred.dat")
                logfile.println("#Prediction on " + refname + " with " + ipd)
                val sequence = readSequenceAsString(fastapath, refname)
                val ipds = readInputAsArray(ipdpath + ipd)
                val cpgs = findAllCpG(sequence)

                println("seq : %d\nipd : %d\ncpg : %d"
                        .format(sequence.length, ipds.length, cpgs.length))

                val prfs: List[(Int, Profile)] = cpgs.map(i => (i, pointProf(ipds, i))).filter(_._2.mcv > 0)
                val ita = prfstoIta(prfs) // use your beta (task parameter ?)
                
                println("len of ita = %d".format(ita.size))

                val optseg = callSegmentationIta(prfs, gamma, ita)
                val scores = if (optseg.nonEmpty) {
                  val ss = segmentSegments(optseg)
                  scoreSegments(ss, ((a:Double) => a))
                  .zip(ss.map(_.length))
                  .map { x => (x._1._1, x._1._2, x._1._3, x._2)}
                } else { List.empty[(Int, Int, Double, Int)] }

                // write out result
                for ( ((a,b,c),p) <- optseg.zip(prfs)) {
                  // idx, ita, class, coverage per strand
                  logfile.println("%d %f %d %f".format(a, b, c, p._2.mcv))
                }
                logfile.close
                val scrfile = new java.io.PrintWriter(outfile + refname + ".avgscr.dat")
                for ( (a,b,c,s) <- scores) {
                  // begin, end, average score, size
                  scrfile.println("%d %d %f %d".format(a, b, c, s))
                }
                scrfile.close

                IOManager.writeRegionsToGFF(outfile, refname, scores)
                IOManager.writeClassToWig(outfile, refname, optseg)
              }
            case _ => println("invalid input file list: %s".format(input)) 
          }
        }
        }
        case _ => println("missing or invalid input(output) field"); sys.exit
      }
    }
    */

    def parseOpt(opts: List[String]): Map[String, String] = opts match {
      case ("-h" :: rest) => Map("help" -> "True")
      case ("--legacy" :: jsonpath :: rest) => Map("legacy" -> jsonpath)

      case ("-i" :: infile :: rest) => Map("inputFile" -> infile) ++ parseOpt(rest)
      case ("-o" :: outfile :: rest) => Map("outputFile" -> outfile) ++ parseOpt(rest)
      case ("-f" :: fasta :: rest) => Map("fastaFile" -> fasta) ++ parseOpt(rest)

      case ("-w" :: wigdir :: rest) => Map("wigDir" -> wigdir) ++ parseOpt(rest)
      case ("--bed" :: bedfile :: rest) => Map("bedFile" -> bedfile) ++ parseOpt(rest)

      case ("-l" :: min_length :: rest) => Map("min_length" -> min_length) ++ parseOpt(rest)
      case ("-b" :: beta :: rest) => Map("beta" -> beta) ++ parseOpt(rest)
      case ("-g" :: gamma :: rest) => Map("gamma" -> gamma) ++ parseOpt(rest)

      case ("-c" :: rest) => Map("continuous" -> "True") ++ parseOpt(rest)

      // case ("-w" :: infile :: rest) => Map("wigdir" -> infile) ++ parseOpt(rest)
      case (command :: Nil) => Map("command" -> command)
      case _ => {
        error("ERROR: parseOpt: ill-formed arguments")
        Map("help" -> "True")
      }
    }

    // of course, deprecated
    def legacyJob(jsonpath: String): Unit = {
      JSON.parseFull(Source.fromFile("task/%s.json".format(jsonpath)).getLines.foldLeft("")(_+"\n"+_)) match {
        case Some(p) => {
          val protocol = p.asInstanceOf[Map[String, Any]]
          protocol("task") match {
            // TODO: Uncomment if you need this
            /*
            case "extract-profile" => extractProfile(protocol)
            case "train-beta0" => trainBeta(protocol)
            case "pointwise-prediction" => pointwisePrediction(protocol)
            case "testOnUnseen" => testOnUnseen(protocol)
            case "predictOnRealCase" => predictOnRealCase(protocol)
            case "reportCoverage" => reportCoverage(protocol)
            */
            case _ => println("unknown task")
          }
          println("task done")
        }
        case _ => println("json parse error")
      }
    }

    // true entry point
    // TODO: Use Monadic pattern
    val opts = parseOpt(args.toList) + ("Full-Commands" -> args.mkString(" "))
    opts.get("legacy") match {
    case Some(path) => legacyJob(path)

    case None =>  opts.get("help") match {
    case Some(_) => Tasks.printHelp()

    case None => opts.get("command") match {
    case Some(command) => command match {
      case "profile" => Tasks.profile(opts)
      case "predict" => Tasks.predict(opts)
      case "makeROC" => Tasks.makeROC(opts)
      case _ => error("unknown command: %s".format(command))
      }
    case None => error("missing command")
    }}}
  }
}
