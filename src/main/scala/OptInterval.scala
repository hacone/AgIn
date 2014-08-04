package OptInterval

import scala.annotation.tailrec

object OptInterval {
  // arg: data, score functtions, minimum length
  def apply[T](data: List[T], scoreFun: List[T=>Double], minLength: List[Int]): List[Int] = {
    require(scoreFun.length == minLength.length,
            "scoreFun and minLength must have same length")
    if (data.isEmpty) return List.empty[Int]

    val numClass = scoreFun.length
    val cumulative = Array.ofDim[Double](data.length+1, numClass)
    val maxScore = Array.ofDim[Double](data.length+1, numClass)
    val tbTable = Array.ofDim[(Int, Int)](data.length+1, numClass)

    // initialise
    for (c <- 0 until numClass) {
      // idx 0 of arrays are sentinels
      cumulative(0)(c) = 0.0
      maxScore(0)(c) = 0.0
      tbTable(0)(c) = (-1, c)
    }

    // recur
    @tailrec
    def rec(d: List[T], i: Int) {
      for (c <- 0 until numClass) {
        cumulative(i)(c) = cumulative(i-1)(c) + scoreFun(c)(d.head)
        if (i < minLength(c)) {
          maxScore(i)(c) = Double.NegativeInfinity
          tbTable(i)(c) = (-2, 0) // must not refer here
        } else {
          // find max for (i,c)
          val el = maxScore(i-1)(c) + scoreFun(c)(d.head) // score for elongated segment
          val ml = {for (cc <- 0 until numClass) yield { // score added for new min len segment
            (maxScore(i-minLength(c))(cc), cc)
          }}.max
          val ns = cumulative(i)(c) - cumulative(i-minLength(c))(c) // score on new segment
          if (ml._1 + ns > el) {
            maxScore(i)(c) = ml._1 + ns
            tbTable(i)(c) = (i-minLength(c), ml._2)
          } else {
            maxScore(i)(c) = el
            tbTable(i)(c) = (i-1, c)
          }
        }
      }
      if (d.tail == Nil) return else rec(d.tail, i+1)
    }
    rec(data, 1)

    // terminate
    val optSeg = scala.collection.mutable.ListBuffer.empty[Int]
    val lastClass = {for (c <- 0 until numClass) yield {
      (maxScore(data.length)(c), c)
    }}.max._2

    var now = (data.length, lastClass)
    var next = (data.length, lastClass)
    while (next._1 >= 0) {
      if (now._1 == next._1) {
        next._2 +=: optSeg
        now = (next._1 - 1, next._2)
        next = tbTable(next._1)(next._2)
        // println("i:%d c:%d".format(next._1, next._2))
      } else {
        now._2 +=: optSeg
        now = (now._1 - 1, now._2)
      }
    }
    /*
    if (next._1 == -2)
      sys.error("something wrong has happened in segmentation")
    else */
      optSeg.tail.result
  }
}
