package hacone.AgIn
import Classifier._

object PrintData {
    def printProfVsBis(prfs: List[(Int, Profile)], bis: Array[Bisulfite], outfile: String) {
      val vsfile = new java.io.PrintWriter(new java.io.FileWriter(outfile, true))
      var i = 0
      for (s <- prfs) {
        while (bis(i).position < s._1) {
          i += 1
          if (i == bis.length) return ()
        }
        if (bis(i).position == s._1) {
          if (bis(i).score >= 0 && s._2.mcv > 0) {
            var str = "%d,".format(i)
            for (j <- -10 to 10) str += "%f,".format(s._2(j))
            str += {
              if (bis(i).score < 0.5) "1" else "0"
            }
            vsfile.println(str)
          }
        }
      }
      vsfile.close
    }
}
