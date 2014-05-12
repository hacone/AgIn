# make package for distribution
bin/sbt package-dist
sed -ie "s/utgenome.sample.Hello/hacone.AgIn.AgIn/" target/dist/bin/classworld.conf
