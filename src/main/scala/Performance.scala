package hacone.AgIn

// Compare a pointwise character against a regional character

trait Performance {
  def eval(prediction: List[(Int, Int)]): Map[String, Double]
}

class Performance_Bed(val bedfile: String) extends Performance {

  val reg = scala.io.Source.fromFile(bedfile)
            .getLines.toList.map { line =>
    val ls = line.split(",")
    (ls(1).toInt, ls(2).toInt, ls(3).toDouble)
  }

  type Bool = Boolean
  def compareRegion[T,U,V](regions: List[(Int, Int, T)], char: List[(Int, U)],
                           fun: (T,U) => V): List[(Int, V)] = {
    val lb = scala.collection.mutable.ListBuffer.empty[(Int, V)]

    @scala.annotation.tailrec
    def loop(r: List[(Int, Int, T)], c: List[(Int, U)]): Unit = {
      if (r.isEmpty || c.isEmpty) { return }
      else if (c.head._1 < r.head._1) { // missing corresponding region
        loop(r, c.tail) 
      } else if (c.head._1 <= r.head._2) {
        lb += ((c.head._1, fun(r.head._3, c.head._2)))
        loop(r, c.tail)
      } else loop(r.tail, c)
    }
    loop(regions, char)
    lb.result
  }

  sealed abstract case class Prediction(i: Int)
  object TP extends Prediction(0)
  object FP extends Prediction(1)
  object FN extends Prediction(2)
  object TN extends Prediction(3)
  object Pass extends Prediction(4)

  // class 0 means Hypo-methylation == Positive
  def ff(x: Double, y: Int): Prediction = {
    if (x < 0.3 && y == 0) TP
    else if (x < 0.3 && y == 1) FN
    else if (x > 0.7 && y == 0) FP
    else if (x > 0.7 && y == 1) TN
    else Pass
  }

  // count up TP, FN, FP, TNs
  def eval(prediction: List[(Int, Int)]): Map[String, Double] = {

    // TODO: fusion
    val result = compareRegion(reg, prediction, ff)
    val tps = result.count(_._2 == TP)
    val fns = result.count(_._2 == FN)
    val fps = result.count(_._2 == FP)
    val tns = result.count(_._2 == TN)

    val mccDenom = math.sqrt((tps+fps).toDouble) * math.sqrt((tps+fns).toDouble) * 
                   math.sqrt((tns+fps).toDouble) * math.sqrt((tns+fns).toDouble)
    
    Map(
      "TP" -> tps,
      "TN" -> tns,
      "FP" -> fps,
      "FN" -> fns,
      "Sens" -> tps.toDouble / (tps + fns),
      "Spec" -> tns.toDouble / (tns + fps),
      "Prec" -> tps.toDouble / (tps + fps),
      "NPV" -> tns.toDouble / (tns + fns),
      "F1" -> 2 * (tps.toDouble) / (2*tps + fps + fns),
      "MCC" -> ((tps*tns) - (fps*fns)).toDouble / mccDenom,
      "Acc" -> (tns + tps).toDouble / (tps + tns + fps + fns)
      )
    }
}

// TODO: for wig file
class Performance_Wig(val score_wigfile: String, val cover_wigfile: String) extends Performance {

  val bis = IOManager.readWigAsArray(score_wigfile, cover_wigfile)

  type Bool = Boolean
  def compareRegion[T,U,V](bisul: Array[Bisulfite], char: List[(Int, U)],
                           fun: (Bisulfite,U) => V): List[(Int, V)] = {
    val lb = scala.collection.mutable.ListBuffer.empty[(Int, V)]
    @scala.annotation.tailrec
    def loop(r: List[Bisulfite], c: List[(Int, U)]): Unit = {
      // println("r.head.position: %d  c.head._1: %d".format(r.head.position, c.head._1))
      if (r.isEmpty || c.isEmpty) { return }
      else if (c.head._1 == r.head.position) {
        lb += ((c.head._1, fun(r.head, c.head._2)))
        loop(r.tail, c.tail)
      } else if (c.head._1 < r.head.position) { // missing corresponding bisul data
        loop(r, c.tail) 
      } else loop(r.tail, c)
    }
    loop(bisul.toList, char)
    // println("len(bisul)=%d len(char)=%d".format(bisul.length, char.length))
    // println("successfully compared: %d".format(lb.result.length))
    lb.result
  }

  sealed abstract case class Prediction(i: Int)
  object TP extends Prediction(0)
  object FP extends Prediction(1)
  object FN extends Prediction(2)
  object TN extends Prediction(3)
  object Pass extends Prediction(4)

  // class 0 means Hypo-methylation == Positive
  def ff(x: Bisulfite, y: Int): Prediction = {
    val sc = x.score
    if (sc < 0.3 && y == 0) TP
    else if (sc < 0.3 && y == 1) FN
    else if (sc > 0.7 && y == 0) FP
    else if (sc > 0.7 && y == 1) TN
    else Pass
  }

  // count up TP, FN, FP, TNs
  def eval(prediction: List[(Int, Int)]): Map[String, Double] = {

    // TODO: fusion
    val result = compareRegion(bis, prediction, ff)
    val tps = result.count(_._2 == TP)
    val fns = result.count(_._2 == FN)
    val fps = result.count(_._2 == FP)
    val tns = result.count(_._2 == TN)

    // println("TP:%d FN:%d FP:%d TN:%d".format(tps, fns, fps, tns))

    val mccDenom = math.sqrt((tps+fps).toDouble) * math.sqrt((tps+fns).toDouble) * 
                   math.sqrt((tns+fps).toDouble) * math.sqrt((tns+fns).toDouble)
    
    Map(
      "TP" -> tps,
      "TN" -> tns,
      "FP" -> fps,
      "FN" -> fns,
      "Sens" -> tps.toDouble / (tps + fns),
      "Spec" -> tns.toDouble / (tns + fps),
      "Prec" -> tps.toDouble / (tps + fps),
      "NPV" -> tns.toDouble / (tns + fns),
      "F1" -> 2 * (tps.toDouble) / (2*tps + fps + fns),
      "MCC" -> ((tps*tns) - (fps*fns)).toDouble / mccDenom,
      "Acc" -> (tns + tps).toDouble / (tps + tns + fps + fns)
    )
  }
}
