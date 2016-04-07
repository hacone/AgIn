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

      case (command :: Nil) => Map("command" -> command)
      case _ => {
        error("ERROR: parseOpt: ill-formed arguments")
        Map("help" -> "True")
      }
    }

    // true entry point
    // TODO: Use Monadic pattern
    val opts = parseOpt(args.toList) + ("Full-Commands" -> args.mkString(" "))

    opts.get("help") match {
    case Some(_) => Tasks.printHelp()

    case None => opts.get("command") match {
    case Some(command) => command match {
      case "profile" => Tasks.profile(opts)
      case "predict" => Tasks.predict(opts)
      case "makeROC" => Tasks.makeROC(opts)
      case "showScores" => Tasks.showScores(opts)
      case "profileOn2mer" => Tasks.profileOn2mer(opts)
      case _ => error("unknown command: %s".format(command))
      }
    case None => error("missing command")
    }}
  }
}
