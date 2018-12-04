package utgenome

import sbt._
import Keys._

object ProjectBuild extends Build {
   lazy val root = Project(
     id ="AgIn",  // Set your project name here
     base = file("."),
     settings = 
       Defaults.defaultSettings 
       ++ Seq(PackageTask.packageDistTask) 
       ++ PackageTask.distSettings 
       ++ Seq(
       	  scalaVersion := "2.10.7",
          organization := "ysuzuki.mlab",
       	  version := "1.0",
       	  scalacOptions ++= Seq("-encoding", "UTF-8", "-deprecation", "-unchecked"),
          javaOptions += "-Xmx8G",
          parallelExecution := true,
          crossPaths := false,
       	  libraryDependencies ++= Seq(
            "org.codehaus.plexus" % "plexus-classworlds" % "2.4",
            "org.scalatest" % "scalatest_2.10" % "2.2.0" % "test",
	          // Add other libraries here
            "net.sf.opencsv" % "opencsv" % "2.3",
            // "io.argonaut" %% "argonaut" % "6.1-M2",
            // "io.argonaut" %% "argonaut" % "6.0.4",
            // "org.scalaz" % "scalaz-core_2.11" % "7.1.0-M7",
            "org.xerial" % "xerial-core" % "3.2.2"
	        )
        )
   )
}
