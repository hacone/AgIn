package hacone.AgIn.GridSearcher

object GridSearcher {
    def apply(l: (Double, Double), s: Int): List[Double] = {
        // require(s > 0 && l._1 < l._2, "empty grid")
        val w = ((l._2 - l._1) / s)
        val ll = { for (i <- 0 to s) yield l._1 + w * i }
        ll.toList
    }
    def next(l: List[(Double, Double)], d: Double): List[Double] = {
        require(l.nonEmpty, "empty grid")
        val s = l.length
        val ls = l.min._1
        val le = l.max._1
        val w = (le - ls) / s.toDouble

        if (l.foldLeft(true)((t,v)=>t && v._2==l.head._2)) {
            println("expanding search space")
            return {for (i <- 0 to s) yield (ls*1.5 - le*0.5) + w*i*2}.toList
        }

        val (ms, mi) = {for ((a,b) <- l; if !b.isNaN) yield (b,a)}.max
        val (ns, ni) = {for ((a,b) <- l; if !b.isNaN) yield (b,a)}.min

        println("mi = %f".format(mi))
        if (ms - ns < d) return List(mi)

        if (mi == le) {
            return {for (i <- 0 to s) yield (ls*1.5 - le*0.5) + w*i}.toList
        } else if (mi == ls) {
            return {for (i <- 0 to s) yield (ls*0.5 + le*0.5) + w*i}.toList
        } else {
            // return {for (i <- 0 to s) yield mi - (le-ls)*0.125 + w*i*0.25}.toList
            return {for (i <- 0 to s) yield mi - (le-ls)*0.25 + w*i*0.5}.toList
        }
    }
}
