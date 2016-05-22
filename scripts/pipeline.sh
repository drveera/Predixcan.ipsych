#!/bin/sh

OPTS=`getopt -o h -l predict,extract,assoc,manhattan,help,outdir:,chunkdir:,pheno:,weights:,\
results:,cov:,expression:,infoscore:,linear:,dosagefolder:,alldbs,dbsdir:  -n "predixcan" -- "$@"`

eval set -- "$OPTS"
helpmessage() {
echo "For extraction
usage: predixcan --extract --outdir <value> --chunkdir <value> 

For prediction
usage: predixcan --predict --outdir <value>  --weights <value>  
note: the outdir should contain a folder 'dosages' with extracted dosage files

For association
usage: predixcan --assoc --outdir <value>  --pheno <value> --cov <value>
note: the outdir should contain a folder 'expression' with predicted expressions 

For manhattan plots
usage: predixcan --manhattan --outdir <value> 
Note: the outdir should contain a folder 'results' with association results. 

First three can be done at once:
usage: predixcan --extract --predict --assoc  --outdir <value> --chunkdir <value> --weights <value> --pheno <value>  

Manhattan plotting should be done seperately since they are not submitted as jobs to cluster

arguments details:

logical arguments that takes **no values**:
--extract 
--predict 
--assoc 
--manhattan 
--linear


arguments that require values:
--outdir <name of the output folder with full path>
--chunkdir <full path of the qc1 folder containing the dosage files>
--dosagefolder <full path of the folder containing extracted dosages>
--weights <full path of the folder containing the precalculated weights (.db files) >
--pheno <phenotype file>
--cov <covariate file> [this is optional, so can be skipped if not required]
--expression <full path of the folder containing the predicted expressions>
--results <full path of the folder containing the association results> 
"

}

while true; do
	case "$1" in
	--predict) predict=true; shift ;;
	--extract) extract=true; shift;;
	--assoc) assoc=true; shift;;
	--manhattan) manhattan=true; shift;;
	--outdir) outdir=$2;  shift 2;;
	--chunkdir) chunkdir=$2; shift 2;;
	--pheno) pheno=$2; shift 2;;
	--pred_exp) pred_exp=$2; shift 2;;
	--weights) weights=$2; shift 2;;
	--dosagefolder) dosedir=$2; shift 2;;
	--expression) expression=$2; shift 2;;
	--results) results=$2; shift 2;;
	--cov) cov_file=$2; shift 2 ;;
	--linear) linear=true; shift ;;
	--alldbs) alldbs=true; shift ;;
	--dbsdir) dbsdir=$2; shift 2;;
	-h | --help) helpmessage; exit 1 ; shift;;
	--) echo "type -h for help"; shift; break;;
	*) break ;;
	esac
done


scriptdir=`dirname $0`
#check if alldbs given, if so then submit full jobs

if [ `echo $alldbs` ];
then	
	check_outdir
	ls -1 $dbsdir/*/ > $outdir/dbslist

	sh $scriptdir/pipeline.sh --predict 
fi


# check if there is atleast one  of the four modes is given

if !([ `echo "$extract"` ] || [ `echo "$predict"` ] || [ `echo "$assoc"` ] || [ `echo "$manhattan"` ]) ;
then
	echo "atleast one of the logical arguments is required; see help "
	exit 1
fi


# evaluation functions

check_outdir (){
	if [ -z $outdir ] ;
	then
		echo "--outdir is required"
		exit 1
	fi
}

check_chunkdir (){
	if [ -z $chunkdir ];
	then
		echo "--chundir is required"
		exit 1
	fi
}

check_dosedir (){
	if [ -z $dosedir ];
	then
		echo "--dosedir is required"
		exit 1
	fi
}

check_weights () {
	if [ -z $weights ];
	then
		echo "--weights is required"
		exit 1
	fi
}

check_pheno () {
	if [ -z $pheno ];
	then
		echo "--pheno is required"
		exit 1
	fi
}

check_weightsfolder () {
	if [ ! -d $weights ];
	then
		echo "weights should be a folder"
		exit 1
	fi
}

#if --dosagefolder is not specified, save the dir inside $dosedir
if [ -z $dosedir ];
then 
	dosedir=$outdir/dosages
fi

check_dosagefolder () {
	if [ ! -d $dosedir ];
	then
		echo "either the 'dosages' folder not found inside $outdir [or] the dosage folder you specified doesn't exist"
		exit 1
	fi 
}

check_dosageinside () {
	ls $dosedir/* > $outdir/dosages.list
	if [ ! -s $outdir/dosages.list ];
	then
		echo "the dosage folder is empty"
		exit 1
	fi
	
}

check_expressionfolder () {
	if [ ! -d $outdir/expression ];
	then
		echo "folder 'expression' not found inside $outdir"
		exit 1
	fi
}

check_resultsfolder () {
	if [ ! -d $outdir/results ];
	then
		echo "folder 'results' not found inside $outdir"
		exit 1
	fi	
}

check_expressionlist () {
	if [ ! -s $outdir/expression.list ];
	then
		echo "folder expression is empty!"
		exit 1
	fi
}
check_resultslist () {
	if [ ! -s $outdir/results.list ];
	then
		echo "folder results is empty"
		exit 1
	fi
}

check_samples() {
	if [ ! -s $dosedir/samples.txt ];
	then
		echo "samples.txt is not found inside dosage folder"
		exit 1
	fi
}

check_uresultsfolder() {
	if [ ! -d $results ];
	then
		echo "the argument --results require a folder (not a file)"
		exit 1
	fi
}

check_uresultsinside(){ #u for user defined
	ls $results/* > $outdir/uresults.list
	if [ ! -s $outdir/uresults.list ] ;
	then
		echo "the folder you specified for --results is empty!"
		exit 1
	fi
	rm $outdir/uresults.list
}

check_uexpressionfolder() {
	if [ ! -d $expression ];
	then
		echo "the argument --expression require a folder(not a file)"
		exit 1
	fi
}

check_uexpressioninside() {
	ls $expression/* > $outdir/uexpression.list
	if [ ! -s $outdir/uexpression.list ];
	then
		echo "the folder you specified for --expression is empty"
		exit 1
	fi

}

makejobid () {
         cat $outdir/job${1}ID | awk '{print $4}'
        }

#echo "exiting"
#exit 1

#script variables
#chunkdir
#outdir
#pheno
#weights
scriptdir=`dirname $0`
softdir=`dirname $scriptdir`

#Required softwares
source /com/extra/R/LATEST/load.sh



#EXTRACTION **********

if [ `echo "$extract"` ];
then
	#check arguments
	check_outdir 
	check_chunkdir
	
	maindir=`dirname $chunkdir`

	#make a list of chunks chromosome wise
	for cr in `seq 22`
	do
		ls $chunkdir/*chr${cr}_*gz > $outdir/chr$cr.chunks.list
		ls $maindir/info/*chr${cr}_*info > $outdir/chr$cr.info.list
	done




	#write job script to extract hapmap snps from all chunks
	for cr in `seq 22`
	do
		while read chunk <&3 && read info <&4 
		do
		subsetname=`basename $chunk`
		echo "Rscript $softdir/scripts/extract.hapmap.R \
			$softdir/resources/hapmap/chr$cr.chr.txt \
			$chunk \
			$info \
			$outdir/$subsetname.subset" 
		done 3<$outdir/chr$cr.chunks.list 4<$outdir/chr$cr.info.list
	done > $outdir/job1.adispatch

	#submit the job
	head -1 $outdir/job1.adispatch > $outdir/testjob1.adispatch
	adispatch --mem=8g $outdir/testjob1.adispatch > $outdir/job1aID
	job1aID=`makejobid 1a`
	adispatch -p normal --mem=8g --dependency=afterok:$job1aID $outdir/job1.adispatch > $outdir/job1ID
	job1ID=`makejobid 1`

	#once the extraction is finished, concatenate dosage chunks chromsosome wise
	mkdir $outdir/dosages

	for i in `seq 22`
	do
		echo "cat $outdir/*chr${i}_*subset > $outdir/dosages/chr$i.dosage"
	done > $outdir/job2.adispatch

	adispatch -p normal --mem=8g --dependency=afterok:$job1ID $outdir/job2.adispatch > $outdir/job2ID
	job2ID=`makejobid 2`

	#submit jobscript to gzip the dosaage files
	for i in `seq 22`
	do
		echo "gzip $outdir/dosages/chr$i.dosage"

	done > $outdir/job3.adispatch
	adispatch -p normal --mem=8g --dependency=afterok:$job2ID $outdir/job3.adispatch > $outdir/job3ID
	job3ID=`makejobid 3`

	#remove the chunk files
	echo "rm $outdir/*subset
	mkdir $outdir/jobfiles
	cat $outdir/job1*out > $outdir/jobfiles/job1.out 
	rm $outdir/job1*out
	cat $outdir/job2*out > $outdir/jobfiles/job2.out
	rm $outdir/job2*out
	cat $outdir/job3*out > $outdir/jobfiles/job3.out
	rm $outdir/job3*out
	rm $outdir/chr*list" > $outdir/job4.adispatch
	adispatch -p normal --dependency=afterok:$job3ID $outdir/job4.adispatch 
	
fi

#PREDICTION **********

if [ `echo "$predict"` ] ;
then
	#check arguments
	check_outdir
	check_weights
	check_weightsfolder
	check_dosagefolder
	check_dosageinside
	check_samples
	#submit jobscript to predict the gene expression
	mkdir $outdir/expression

	predixcall (){
		echo "python $softdir/scripts/PrediXcan.py --predict \
			--weights $1 \
			--dosages $dosedir \
			--samples samples.txt \
			--output_dir $outdir/expression \
			--outname `basename $1`"
	}
	
	ls $weights/* > $outdir/weights.list

	while read weight
	do
		predixcall $weight 
	done < $outdir/weights.list > $outdir/job5.adispatch

	if [ ! -z $job3ID ];
	then
		adispatch -p normal --mem=64g  --dependency=afterok:$job3ID $outdir/job5.adispatch > $outdir/job5ID
		job5ID=`makejobid 5`
	else 
		adispatch -p normal --mem=64g  $outdir/job5.adispatch > $outdir/job5ID
		job5ID=`makejobid 5`
	fi
fi


#ASSOCIATION **************


if [ `echo "$assoc"` ];
then
	#check arguments
	check_outdir
	check_pheno
	check_expressionfolder
	

	mkdir $outdir/results
	#submit the jobscript to run association
	if [ `echo "$linear"` ];
	then
		testtype=linear
	else
		testtype=logistic
	fi
	
	if [ -z $cov_file ];
	then
		assoc_call(){
			echo "python $softdir/scripts/PrediXcan.py --assoc \
			--pheno $pheno \
			--pred_exp $1 \
			--output_dir $outdir/results \
			--outname `basename $1` \
			--$testtype \
			--nthread 8
			"
		}
	else
	        assoc_call(){
                        echo "python $softdir/scripts/PrediXcan.py --assoc \
                        --pheno $pheno \
                        --pred_exp $1 \
                        --output_dir $outdir/results \
                        --outname `basename $1` \
			--cov $cov_file \
			--$testtype \
			--nthread 8 
                        "
                }
	fi
	if [ ! -z $expression ];
	then
		check_uexpressionfolder
		ls $expression/* > $outdir/expression.list
	else
		ls $outdir/expression/* > $outdir/expression.list
	fi
	check_expressionlist
	
	while read i
	do
		assoc_call $i
	done < $outdir/expression.list > $outdir/job6.adispatch
	
	if [ ! -z ${job5ID+x} ];
	then
		adispatch -p normal -c 8 --mem=256g --dependency=afterok:$job5ID $outdir/job6.adispatch > $outdir/job6ID
		job6ID=`makejobid 6`
	else 
		adispatch -p normal -c 8 --mem=256g $outdir/job6.adispatch > $outdir/job6ID
		job6ID=`makejobid 6`
	fi

fi

#MANHATTAN ************

if [ `echo "$manhattan"` ] ;
then
	#check arguments
	check_outdir
	check_resultsfolder
	

	mkdir $outdir/plots
	
	manhattan_call() {
		Rscript $softdir/scripts/manhattan.R $1 $softdir/resources/ensembl.bed $2
	}
	
	if [ ! -z $results ];
	then
		check_uresultsfolder
		ls $results/* > $outdir/results.list
	else
		ls $outdir/results/* > $outdir/results.list
	fi
	check_resultslist
	
	while read i
	outname=`basename $i`
	do
		manhattan_call $i $outdir/plots/$outname
	done < $outdir/results.list
fi


