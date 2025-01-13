#!/bin/bash

# set -xv

# to retrieve insertion sites

Rviz=$1


readslist=`cat $Rviz  | sed 's/,[^\t]*//g' | awk '{OFS="\t"}{if ($15~/60/){print $1,$2,$3,$7,$8,$9,$10,$11,$21";"$23";"$24";"$25";"$30}}' | cut -f1  | sort | uniq`



> $Rviz".insertions.bed"

for i in $readslist; do
    instance=$(mktemp ./instance.XXXXXX)
    cat $Rviz  | sed 's/,[^\t]*//g' | awk '{OFS="\t"}{print $1,$2,$3,$7,$8,$9,$10,$11,$30,$21";"$24";"$25";"$30}' | grep $i  | sed 's/ATCOPIA93.2/ATCOPIA93/g' > $instance
    read=`cat $instance | cut -f1 | uniq`
    arrangement=`cat $instance | cut -f9 | tr "\n" "|"`
    orientation=`cat $instance | cut -f8 | tr "\n" "|"`
    INSA=`cat $instance | grep "\."  | cut -f5  | head -1`
    INSB=`cat $instance | grep "\." | cut -f6  | head -1`
    CHR=`cat $instance | grep "\." | cut -f4  | head -1`
    TEtype=`cat $instance | cut -f10  | grep -v "\.$"`

    echo $i


    #  echo $read
    #  echo ""
    #  echo $arrangement
    #  echo ""
    #  echo $orientation
    #  echo ""
    #  echo $INSA
    #  echo ""
    #  echo $INSB
    #  echo ""
    #  echo $CHR
    #  echo ""
    #  echo $TEtype
    #  exit 1

    ## +|+| = +   ; -|-| = +   ; +|+|+ = +     ; -|-|- = +
    if [ $orientation = "+|+|" ] || [ $orientation = "-|-|" ]  || [ $orientation = "+|+|+" ] || [ $orientation = "-|-|-" ]; then
        INSSTRAND="+"
    ## +|-| = -   ; +|- = -  ; +|-|+ = -   ; -|+|- = -
    elif [ $orientation = "+|-|" ] || [ $orientation = "-|+|" ]  || [ $orientation = "+|-|+" ] || [ $orientation = "-|+|-" ]; then
        INSSTRAND="-"
    fi

    # checking if arrangement is cx
    if [[ "$arrangement" =~ ^\.\|[^\|]+\|\.\|$ ]]; then
        echo "arrangement is cx NOTE - TE - NOTE" 
        if [[ "$orientation" =~ ^\+\|[^\|]+\|\+\|$ ]]; then
        # orientation is +
        # arrangement is cx+
            RESULT=$CHR"\t"${INSB}"\t"${INSB}"\t"${INSSTRAND}"\t"${TEtype}"\t"${read}
        elif [[ "$orientation" =~ ^\-\|[^\|]+\|\-\|$ ]]; then
        # orientation is -
        # arrangement is cx-
            RESULT=$CHR"\t"${INSA}"\t"${INSA}"\t"${INSSTRAND}"\t"${TEtype}"\t"${read}
        else
            echo "${instance} -> something is wrong here1"
            RESULT=

        fi
    
    # checking if arrengement is notete
    elif [[ $arrangement =~ ^\.\|[^\|]+|\$ ]]; then
    #  exit 1
        echo "arrangement is notete"
        if [[ $orientation =~ ^\+\|[^\|]+\|$ ]]; then
        # orientation is +
        # arrangement is notete+
            RESULT=$CHR"\t"${INSB}"\t"${INSB}"\t"${INSSTRAND}"\t"${TEtype}"\t"${read}
        elif [[ $orientation =~ ^\-\|[^\|]+\|$ ]]; then
        # orientation is -
        # arrangement is notete-
            RESULT=$CHR"\t"${INSA}"\t"${INSA}"\t"${INSSTRAND}"\t"${TEtype}"\t"${read}
        else
            echo "${instance} -> something is wrong here2"
            RESULT=
        fi

    # checking if arrengement is tenote
    elif [[ $arrangement =~ ^[^\|^\.]+\|\.\|$ ]]; then
        echo "arrangement is tenote"
        if [[ $orientation =~ ^[^\|]+\|\+\|$ ]]; then
        # orientation is +
        # arrangement is tenote+
            RESULT=$CHR"\t"${INSA}"\t"${INSA}"\t"${INSSTRAND}"\t"${TEtype}"\t"${read}
        elif [[ $orientation =~ ^[^\|]+\|\-\|$  ]]; then
        # orientation is -
        # arrangement is tenote-
            RESULT=$CHR"\t"${INSB}"\t"${INSB}"\t"${INSSTRAND}"\t"${TEtype}"\t"${read}
        else
            echo "${instance} -> something is wrong here3"
            RESULT=
        fi
    else
        echo "${instance} -> something is VERY wrong here"
        RESULT="# DIFFICULT CASE TO AUTOMATE: "`cat $Rviz  | sed 's/,[^\t]*//g' | awk '{OFS="\t"}{print $1,$2,$3,$7,$8,$9,$10,$11,$30}' | grep $i`
    
    

    fi
    echo -e $RESULT | cut -f 1-6 >> $Rviz".insertions.bed"
    rm $instance
done

sed -i '/^$/d' $Rviz".insertions.bed"



grep "^#" $Rviz".insertions.bed" > $Rviz".insertions.cov.bed"
sort -k1,1 -k2,2n $Rviz".insertions.bed" | bedtools merge -i stdin -c 4,5,5,6 -o collapse,count,collapse,collapse >> $Rviz".insertions.cov.bed"