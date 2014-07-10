
package Classifier {

  // profile: List[(relative pos, mean ipd, cov, var, #pos)]
  case class Profile(val profile : List[(Int, Double, Int, Double, Int)]) {

    val (prof, meancov) = {
      val b = new Array[(Double, Int, Double, Int)](21)
      var s = 0
      var n = 0
      for (p <- profile) {
        b(p._1 + 10) = (p._2, p._3, p._4, p._5)
        s += p._3 // p._3 is sum of coverage = (fc + rc)
        n += p._5 // n = 2
      }
      (b, if (n > 0) s/(n.toDouble) else 0.0)
    }

    def apply(i: Int): Double = prof(i+10)._1
    def cov(i: Int): Int = prof(i+10)._2
    def pos(i: Int): Int = prof(i+10)._4
    def mcv(): Double = meancov

    // def toString: String = 
  }

  object Predictor {
    def makeIta(prfs: List[(Int, Profile)], vn: List[Double]): List[Double] = {
      def weightedProd(v: List[Double], u: List[Double], w: List[Double]): Double = {
        v.zip(u).zip(w).map{ x: ((Double, Double), Double) =>
          x match { case ((a, b), c) => a*b*c }
        }.sum
      }
      prfs.map { p =>
        weightedProd(
          (-10 to 10).toList.map { i => p._2(i) }, vn,
          (-10 to 10).toList.map { i => p._2.cov(i).toDouble }
        )
      }
    }
  }
}
