#
#  Figure out what project this script belongs to and include the right files
#  Note that variable $project is set for use by the rest of the script
#
project=topmed
if [[ $0 =~ (inpsy) ]]; then
  project=inpsyght
fi
echo "This is project $project"
source /usr/cluster/$project/bin/${project}_defs.inc
source /usr/cluster/$project/bin/${project}_functions.inc

