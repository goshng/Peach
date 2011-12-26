sudo rm -rf /Applications/rnaseq.app
sudo rm -rf /Applications/rnaseq64.app
sudo rm -rf /Library/Frameworks/RNAseq.framework
exit
sudo ./pkg R-2.14.0
exit

sudo rm -rf /Applications/rnaseq.app
sudo rm -rf /Applications/rnaseq64.app
sudo rm -rf /Library/Frameworks/RNAseq.framework
sudo cp -r /Users/goshng/Documents/Projects/peach/rnaseq/build/Development/rnaseq.app /Applications
sudo cp -r /Users/goshng/Documents/Projects/peach/rnaseq/build/Development/rnaseq.app /Applications/rnaseq64.app
sudo cp -r /Library/Frameworks/RNAseq.framework-original /Library/Frameworks/RNAseq.framework
exit

./buildR R-2.14.0 i386 snowleopard-i386
./install-packages.sh
./install-packages.sh
sudo ./universal
sudo rm -rf /Applications/rnaseq.app
sudo rm -rf /Applications/rnaseq64.app
sudo cp -r /Users/goshng/Documents/Projects/peach/rnaseq/build/Development/rnaseq.app /Applications
sudo cp -r /Users/goshng/Documents/Projects/peach/rnaseq/build/Development/rnaseq.app /Applications/rnaseq64.app
chown -Rh root:admin /Applications/rnaseq.app /Applications/rnaseq64.app
chmod -R g+w /Applications/rnaseq.app /Applications/rnaseq64.app

sudo mkdir /Applications/rnaseq.app/Others
sudo cp /Users/goshng/Documents/Projects/peach/rdev/downloads/cutadapt-1.0.tar.gz /Applications/rnaseq.app/Others
sudo cp /Users/goshng/Documents/Projects/peach/rdev/package/bwa /Applications/rnaseq.app/Others
sudo cp /Users/goshng/Documents/Projects/peach/rdev/package/samtools /Applications/rnaseq.app/Others
chown -Rh root:admin /Applications/rnaseq.app /Applications/rnaseq64.app
chmod -R g+w /Applications/rnaseq.app /Applications/rnaseq64.app

sudo ./pkg R-2.14.0
sudo rm -rf /Applications/rnaseq.app
sudo rm -rf /Applications/rnaseq64.app
sudo rm -rf /Library/Frameworks/RNAseq.framework
sudo mv /Library/Frameworks/RNAseq.framework /Library/Frameworks/RNAseq.framework-original
exit


cp /Library/Frameworks/R.framework/Resources/lib/libRblas.vecLib.dylib .
cp /Library/Frameworks/R.framework/Resources/lib/libreadline.5.2.dylib .
cp /Library/Frameworks/R.framework/Resources/lib/libgfortran.2.dylib .
cp /Library/Frameworks/R.framework/Resources/lib/libgcc_s.1.dylib .


cp leopard-x86_64/R-2.14.0-leopard-x86_64-bld.tar.gz deploy/leopard/R-2.14.0/x86_64/R-2.14.0-leopard-universal.tar.gz
cp /Users/goshng/Documents/Projects/peach/rnaseq/build/Development/rnaseq-1.0-Leopard.tar.gz /Builds/R-builds/nightly/deploy/leopard/R-2.14.0/x86_64
cp /Users/goshng/Documents/Projects/peach/rnaseq/build/Development/rnaseq-1.0-Leopard64.tar.gz /Builds/R-builds/nightly/deploy/leopard/R-2.14.0/x86_64
