package hacone.AgIn

import hacone.AgIn._

object Tasks extends xerial.core.log.Logger {
    // TODO: this won't run for now; write this today !
    def profile(opts: Map[String, String]): Unit = {

      // should raise some warning if filenames are not specified
      val infile = opts.get("inputFile").getOrElse("modifications.csv")
      val outfile = opts.get("outputFile").getOrElse("profiles.dat")
      val fasta = opts.get("fastaFile").getOrElse("default.fasta")
      val wigdir = opts.get("wigDir").getOrElse("input/wig/")

      val m_outpw = new java.io.PrintWriter(new java.io.FileWriter("methyl_" + outfile, true))
      val um_outpw = new java.io.PrintWriter(new java.io.FileWriter("unmethyl_" + outfile, true))

      val avg_mprofs = new Array[Double](21)
      val avg_uprofs = new Array[Double](21)
      var n_mprofs = 0
      var n_uprofs = 0

      for ((refname, ipds) <- IOManager.loadIPD(infile)) {
        // info("Handling %s : len(ipds) = %d".format(refname, ipds.length))
        info("Handling %s : len(ipds) = ???".format(refname))

        val dna_seq = IOManager.readSequenceAsString(fasta, refname)
        info("len(dna_seq) = %d".format(dna_seq.length))

        val bisul = IOManager.readWigAsArray(wigdir+refname+".scr.wig", wigdir+refname+".cov.wig")
        info("len(bisul) = %d".format(bisul.length))

        val cpgs = Profiling.findAllCpG(dna_seq)
        info("len(cpgs) = %d".format(cpgs.length))
        info(cpgs.take(10).mkString(","))

        val (mcpgs, ucpgs) = Profiling.partitionMethylationState(cpgs, bisul)
        info("(#M, #Um) = (%d, %d)".format(mcpgs.length, ucpgs.length))

        // TODO: anyway, I can treat dna_seq and cpgs as iterator
        // TODO: how to abstract them ?
        @scala.annotation.tailrec
        def loop(midx: List[Int], uidx: List[Int]): Unit = {

          if (midx.isEmpty && uidx.isEmpty) {
            ()
          } else if (uidx.isEmpty) {
            for (prf <- Profiling.makeProfile(ipds, midx.head); if prf.mcv > 0) {
              m_outpw.println( (-10 to 10).toList.map {i => "%.4f".format(prf(i))}.mkString(" "))
              (-10 to 10).foreach { i => avg_mprofs(i+10) += prf(i) }
              n_mprofs += 1
            }
            loop(midx.tail, uidx)
          } else if (midx.isEmpty) {
            for (prf <- Profiling.makeProfile(ipds, uidx.head); if prf.mcv > 0) {
              um_outpw.println( (-10 to 10).toList.map {i => "%.4f".format(prf(i))}.mkString(" "))
              (-10 to 10).foreach { i => avg_uprofs(i+10) += prf(i) }
              n_uprofs += 1
            }
            loop(midx, uidx.tail)
          } else if (midx.head < uidx.head) {
            for (prf <- Profiling.makeProfile(ipds, midx.head); if prf.mcv > 0) {
              m_outpw.println( (-10 to 10).toList.map {i => "%.4f".format(prf(i))}.mkString(" "))
              (-10 to 10).foreach { i => avg_mprofs(i+10) += prf(i) }
              n_mprofs += 1
            }
            loop(midx.tail, uidx)
          } else {
            for (prf <- Profiling.makeProfile(ipds, uidx.head); if prf.mcv > 0) {
              um_outpw.println( (-10 to 10).toList.map {i => "%.4f".format(prf(i))}.mkString(" "))
              (-10 to 10).foreach { i => avg_uprofs(i+10) += prf(i) }
              n_uprofs += 1
            }
            loop(midx, uidx.tail)
          }
        }

        loop (mcpgs, ucpgs)

        // TODO: compute LDA vector here using matrix calculation
        // TODO: for now, I may use simple average as beta vector
        info(avg_mprofs.map(_ / n_mprofs).mkString(", "))
        info(avg_uprofs.map(_ / n_uprofs).mkString(", "))
      }
      // avg for beta vector
      info(avg_mprofs.map(_ / n_mprofs).mkString(", "))
      info(avg_uprofs.map(_ / n_uprofs).mkString(", "))
      m_outpw.close
      um_outpw.close
    }

    def loadBeta(filename: String): List[Double] = filename match {
      // These 2 for compatibility with older spec
      case "P5C3_HdrR" => resources.Resources.p5c3
      case "LDAVector" => resources.Resources.p4c2

      // This is the main case: please specify beta by file
      case fname => {
        try {
          info("loading Beta vector from file: %s".format(fname))
          val _beta: List[Double] = scala.io.Source.fromFile(fname).getLines.toList.map(_.toDouble)
          assert(_beta.length == 21)

          val norm = scala.math.sqrt(_beta.fold(0.0){ (a:Double,b:Double) => a + b*b })
          _beta.map(_ / norm)
        } catch {
          case (ex: Exception) => {
            info(ex)
            info("Default to LDA vector from P4C2.")
            resources.Resources.p4c2
          }
        }
      }
    }
    
    //  this correspond to <trainbeta> task
    def makeROC(opts: Map[String, String]): Unit = {
      val infile = opts.get("inputFile").getOrElse("modifications.csv")
      val outfile = opts.get("outputFile").getOrElse("performance")
      val fasta = opts.get("fastaFile").getOrElse("default.fasta")
      val min_length = opts.get("min_length").getOrElse("50").toInt
      val beta = loadBeta(opts.get("beta").getOrElse("LDAVector"))

      val (refname, ipds) = IOManager.loadIPD(infile).next()
      info("Handling %s : len(ipds) = ???".format(refname))
      val dna_seq = IOManager.readSequenceAsString(fasta, refname)
      info("len(dna_seq) = %d".format(dna_seq.length))
      val cpgs = Profiling.findAllCpG(dna_seq)
      info("len(cpgs) = %d".format(cpgs.length))

      val perf = opts.get("wigDir") match {
        case Some(wigdir) => {
          info("prediction to be compared with "+wigdir+refname+".scr.wig(.cov.wig)")
          new Performance_Wig(wigdir+refname+".scr.wig", wigdir+refname+".cov.wig")
        }
        case None => {
          info("prediction to be compared with "+opts.get("bedFile").getOrElse("default.bed"))
          new Performance_Bed(opts.get("bedFile").getOrElse("default.bed"))
        }
      }

      // position, ita, coverage
      val ita: List[(Int, Double, Double)] = {for {
        i <- cpgs
        prf <- Profiling.makeProfile(ipds, i)
        if prf.mcv > 0
      } yield (i, Classifier.Predictor.toIta(prf, beta), prf.mcv)}.toList

      info("len(ita) = %d".format(ita.length))

      // TODO

      val pw = new java.io.PrintWriter(outfile)
      pw.println("Gamma TP FP FN TN Sens Spec Prec F1 MCC Acc")

      def icopt(iclist: List[Double]): List[(Double, Double)] = {
        iclist.map { ic =>
          val optseg = callSegmentationIta(ita, ic, min_length)
          val toeval = {for (o <- optseg) yield (o._1, o._3)}.toList
          val evp = perf.eval(toeval)
          pw.println("%.3f %d %d %d %d %.4f %.4f %.4f %.4f %.4f %.4f".format(
            ic,
            evp("TP").toInt, evp("FP").toInt, evp("FN").toInt, evp("TN").toInt,
            evp("Sens"), evp("Spec"), evp("Prec"), evp("F1"),
            evp("MCC"), evp("Acc"))) //  toeval.length.toDouble / cpgs.length))
          (ic, evp("MCC"))
        }
      }

      import hacone.AgIn.GridSearcher._
      var icl = GridSearcher((-8, 5),20)
      val deltai = 0.01
      def mid(a: Double ,b: Double) = {(a+b) * 0.5}

      while (icl.max - icl.min > deltai) {
        println("\n icl search : %.5f to %.5f".format(icl.min, icl.max))
        icl = GridSearcher.next(icopt(icl), 0.003)
        println("optimizing (ic = %f)".format(mid(icl.max,icl.min)))
      }

      pw.close
    }

    def predict(opts: Map[String, String]): Unit = {
      val infile = opts.get("inputFile").getOrElse("modifications.csv")
      val outfile = opts.get("outputFile").getOrElse("predict")
      val fasta = opts.get("fastaFile").getOrElse("default.fasta")
      val min_length = opts.get("min_length").getOrElse("50").toInt
      val gamma = opts.get("gamma").getOrElse("-1.80").toDouble
      val beta = loadBeta(opts.get("beta").getOrElse("LDAVector"))

      // whether to perform continuous prediction using arrays of gammas
      val continuous = opts.get("continuous").getOrElse("False").toBoolean
      val commands = opts.get("Full-Commands").getOrElse("???")

      info("minlen, gamma = %d %f".format(min_length, gamma))

      IOManager.writeHeaderToGFF(outfile, commands)

      for ((refname, ipds) <- IOManager.loadIPD(infile)) {
        // info("Handling %s : len(ipds) = %d".format(refname, ipds.length))
        info("Handling %s : len(ipds) = ???".format(refname))

        // TODO: how can I manage this inter dependency of variable
        //       to get better GC and calculation policies ?
        // TODO: fusion 
        val dna_seq = IOManager.readSequenceAsString(fasta, refname)
        info("len(dna_seq) = %d".format(dna_seq.length))
        val cpgs = Profiling.findAllCpG(dna_seq)
        info("len(cpgs) = %d".format(cpgs.length))

        // position, ita, coverage
        val ita: List[(Int, Double, Double)] = {for {
          i <- cpgs
          prf <- Profiling.makeProfile(ipds, i)
          if prf.mcv > 0
        } yield (i, Classifier.Predictor.toIta(prf, beta), prf.mcv)}.toList

        info("len(ita) = %d".format(ita.length))

        if (ita.length > 0) {
          val optseg = callSegmentationIta(ita, gamma, min_length)
          if (optseg.length > 0) {
            val scores = scoreSegments(optseg)
            // write out result
            IOManager.writeRegionsToGFF(outfile, refname, scores)
            IOManager.writeCoverageToWig(outfile, refname, ita)
            IOManager.writeClassToWig(outfile, refname, optseg)
          }

          // TODO: any rational selection of derivation of the gamma list?
          if (continuous) {
            val perturbations = (6 to -3 by -1).toList.map(_/25.0)
            val gamma_list = perturbations.map(p => gamma*(1+p))
            val optseg_continuous = gamma_list.map {
                g => callSegmentationIta(ita, g, min_length).map(_._3)
              }.transpose.map { pos => pos.count(_==1)/10.0 }

            IOManager.writeContinuousPredictionToWig(
              outfile, refname, ita.map(_._1).zip(optseg_continuous))
          }
        }
      }
    }

    def showScores(opts: Map[String, String]): Unit = {
      val infile = opts.get("inputFile").getOrElse("modifications.csv")
      val outfile = opts.get("outputFile").getOrElse("predict")
      val fasta = opts.get("fastaFile").getOrElse("default.fasta")
      // don't need this // val min_length = opts.get("min_length").getOrElse("50").toInt
      val gamma = opts.get("gamma").getOrElse("-1.80").toDouble
      val beta = loadBeta(opts.get("beta").getOrElse("LDAVector"))

      // whether to perform continuous prediction using arrays of gammas
      val continuous = opts.get("continuous").getOrElse("False").toBoolean
      val commands = opts.get("Full-Commands").getOrElse("???")

      info(s"gamma = ${gamma}")

      for ((refname, ipds) <- IOManager.loadIPD(infile)) {
        // info("Handling %s : len(ipds) = %d".format(refname, ipds.length))
        info("Handling %s : len(ipds) = ???".format(refname))

        // TODO: how can I manage this inter dependency of variable
        //       to get better GC and calculation policies ?
        // TODO: fusion 
        val dna_seq = IOManager.readSequenceAsString(fasta, refname)
        info("len(dna_seq) = %d".format(dna_seq.length))
        val cpgs = Profiling.findAllCpG(dna_seq)
        info("len(cpgs) = %d".format(cpgs.length))

        // position, ita, coverage
        val ita: List[(Int, Double, Double)] = {for {
          i <- cpgs
          prf <- Profiling.makeProfile(ipds, i)
          if prf.mcv > 0
        } yield (i, Classifier.Predictor.toIta0(prf, beta), prf.mcv)}.toList

        info("len(ita) = %d".format(ita.length))
        info("# position, ita(score)_raw, coverage")

        // write out in wig

        println(s"chrom=${refname} variableStep")
        ita.foreach { i =>
          println(s"${i._1} ${i._2}")
        }
      }
    }


    def printHelp() {
        println("AgIn: <explanation>")
    }

    // (position, score, class)
    type IDI = (Int, Double, Int)

    // using ita, calculate final scores for given intercept. then call segmentation algorithm
    def callSegmentationIta(ita: List[(Int, Double, Double)], ic: Double, minlen: Int): List[IDI] = {
      // TODO: make it configurable
      val seglen = minlen
      val fm: Double => Double = { s => (s) }
      val fum: Double => Double = { s => -(s) }
      val score = {for ((_,i,c) <- ita) yield i + 2 * c * ic }.toList

      // 0 for unmethylated, 1 for methylated
      val optipdseg = OptInterval.OptInterval(score, List(fum, fm), List(seglen,seglen))

      return {for ((((a, _, _), b), c) <- ita.zip(score).zip(optipdseg)) yield {
        (a, b, c) // position, score, class
      }}.toList
    }

    // TODO: may be good to write this as fold
    // group the elements of a list, from its head, according to result of fn
    def groupWith[T](fn: T => Any)(list: List[T]): List[List[T]] = {
      val b = scala.collection.mutable.ListBuffer.empty[List[T]]

      @scala.annotation.tailrec
      def recur(l: List[T]): Unit = {
        val now = fn(l.head)
        val (hd, tl) = l.span(fn(_) == now)
        b += hd
        if (!tl.isEmpty) recur(tl)
      }

      recur(list)
      b.result
    }


    // get list of (start_idx, end_idx, average score, size) of each predicted region
    def scoreSegments(segs: List[IDI]): List[(Int, Int, Double, Int)] = {
      groupWith((x: IDI) => x._3)(segs).map { seg =>
        val score = seg.map { (x: IDI) => x._2 }.sum / seg.length
        (seg.head._1, seg.last._1, score, seg.size)
      }
    }
}
