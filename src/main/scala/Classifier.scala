import hacone.AgIn.Environment
import scala.collection.immutable.IntMap

package Classifier {
  abstract class Character() { // container for character value
  }
  object Character {
  }

  abstract class Classifier[T <: Character](val f : T => Boolean, val n : String) {
    def apply(c : T) = f(c) 
  }

  case class MeanIPD(val meanipd : Double) extends Character {
  }

  class MeanIPDClassifier(override val f: MeanIPD => Boolean, n : String) extends Classifier[MeanIPD](f, n) {
  }

  // profile: List[(relative pos, mean ipd, cov, var, #pos)]
  case class Profile(val profile : List[(Int, Double, Int, Double, Int)]) extends Character {
    val (profmap, meancov) = {
      val b = Map.newBuilder[Int, (Double, Int, Double, Int)]
      var s = 0
      var n = 0
      for (p <- profile) {
        b += p._1 -> (p._2, p._3, p._4, p._5)
        s += p._3 // p._3 is sum of coverage = (fc + rc)
        n += p._5 // n = 2
      }
      (b.result, if (n > 0) s/(n.toDouble) else 0)
    }

    def apply(i: Int): Double = { profmap(i)._1 }
    def cov(i: Int): Int = { profmap(i)._2 }
    def pos(i: Int): Int = { profmap(i)._4 } // ???
    def mcv(): Double = meancov
  }
  class ProfileClassifier(override val f: Profile => Boolean, n : String) extends Classifier[Profile](f, n) {
  }

  // data: list(character, isPositive)
  object ROC {
    def valueClassifiers[T <: Character](data: List[(T, Boolean)], classifiers: List[Classifier[T]]): List[(Classifier[T], Double, Double, Double, Double, Int, Int)] = {
      for (c <- classifiers) yield {
        var tp, fp, tn, fn = 0
        for (d <- data) {
          (c(d._1), d._2) match {
            case (true, true) => tp += 1
            case (true, false) => fp += 1
            case (false, false) => tn += 1
            case (false, true) => fn += 1
          }
        }
        (c, 
         Sensitivity(tp, fp, tn, fn),
         Specificity(tp, fp, tn, fn),
         PPV(tp, fp, tn, fn),
         NPV(tp, fp, tn, fn), tp, tn
        )
      }
    }

    def Sensitivity(tp : Int, fp : Int, tn : Int, fn : Int) : Double = {
     tp / (tp + fn).toDouble
    }
    def Specificity(tp : Int, fp : Int, tn : Int, fn : Int) : Double = {
     tn / (tn + fp).toDouble
    }
    def PPV(tp : Int, fp : Int, tn : Int, fn : Int) : Double = {
     tp / (tp + fp).toDouble
    }
    def NPV(tp : Int, fp : Int, tn : Int, fn : Int) : Double = {
     tn / (tn + fn).toDouble
    }
  }
}
