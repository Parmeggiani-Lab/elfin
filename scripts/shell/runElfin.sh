: ${execFile='./src/GA/bin/elfin'}
: ${configFile='elfinConfig.json'}
: ${outputDir='output'}
: ${extraArgs=''}

echo Using execFile: $execFile
echo Using configFile: $configFile
echo Using outputDir: $outputDir

cd $HOME/src/elfin

mkdir -p $outputDir
$execFile -s $configFile -o $outputDir $extraArgs
