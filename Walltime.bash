# A function to include in bash scripts to detail walltime for loops


# Get start time
SCRIPTSTART=`date +%s`

# Calling function requires a start time variable (eg walltime $STARTTIME), which is read via the $1 variable
walltime () {
	RED='\033[1;31m'
	GREEN='\033[1;32m'
	YELLOW='\033[1;33m'
	BLUE='\033[1;34m'
	PURPLE='\033[1;35m'
	NOCOLOUR='\033[0m'	
	TIMENOW=`date +%s`
	STARTTIME=$1
	TOTAL=$((TIMENOW-STARTTIME))
	
	if [ "$TOTAL" -lt "120" ]; then
		echo -e "${BLUE}Total wall time: ${PURPLE}$TOTAL ${BLUE}seconds${NOCOLOUR}"
	elif [ "$TOTAL" -lt "7200" ]; then
		TOTAL=$(echo "$TOTAL/60" | bc)
		echo -e "${BLUE}Total wall time: ${PURPLE}$TOTAL ${BLUE}minutes${NOCOLOUR}"
	else
		TOTAL=$(echo "$TOTAL/3600" | bc)
		echo -e "${BLUE}Total wall time: ${PURPLE}$TOTAL ${BLUE}hours${NOCOLOUR}"
	fi
}
