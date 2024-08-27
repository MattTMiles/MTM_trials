#!/bin/bash

#This gets the data and prepares it for Normal and MTM comparison
#The data still needs to be trimmed after for both Normal and MTM so that it shows basically the same thing



psr=$1
#for psr in $(ls -d J*); do 

if [[ ! -f "${psr}/transfer_finished" ]]; then

    echo $psr 
    for dir in $(ls -d /fred/oz005/timing_processed/$psr/20*); do 
        if [[ -f $(echo $dir/*/decimated/*zap*928ch_4p_1t.ar) ]]; then 
            cp $dir/*/decimated/*zap*928ch_4p_1t.ar $psr/data/; 
        elif [[ -f $(echo $dir/*/decimated/*zap.*928ch_4p_1t.ar) ]]; then
            cp $dir/*/decimated/*zap*928ch_4p_1t.ar $psr/data/;
        else
            echo "${dir} doesn't contain file"
            echo "${dir}" >> ${psr}/${psr}_missing.list
        fi; 
    done; 
    #done

    cd $psr/data
    pam -D -e ar_D J*ar

    psradd -T J*ar_D -r 1283.58203125 -o ${psr}_grand_T.ar_D
    cp ${psr}_grand_T.ar_D ..

    cd ../for_timing
    cp ../data/J*ar .

    pam --setnchn 32 -e ar_32ch J*ar

    cd ..

    pam --setnchn 32 -e ar_D_32ch ${psr}_grand_T.ar_D
    psrsmooth -W ${psr}_grand_T.ar_D_32ch

    #Normal timing
    pat -jp -A FDM -C "chan rcvr snr length subint" -f "tempo2 IPTA" -P -s ${psr}_grand_T.ar_D_32ch.sm for_timing/*32ch >> ${psr}.tim_NORMAL

    #Kill bad high SNR obs
    awk '(NR==1) || (NR==2) || ($29 > 15 ) ' ${psr}.tim_NORMAL > temp.tim_NORMAL
    mv temp.tim_NORMAL ${psr}.tim_NORMAL

    #Kill UHF obs
    awk '(NR==1) || (NR==2) || ($13 > 17 ) ' ${psr}.tim_NORMAL > temp.tim_NORMAL
    mv temp.tim_NORMAL ${psr}.tim_NORMAL

    #MTM timing based on Normal timing
    pat -A FDM -p -C "chan rcvr snr length subint" -f "tempo2 IPTA" -P -s /fred/oz002/users/mmiles/MTM_trials/${psr}/${psr}_grand_T.ar_D_32ch.sm $(awk '{print $1}' ${psr}.tim_NORMAL | uniq | cut -d "/" -f2 | grep J | rev | cut -f 2- -d '.' | rev | awk -v awkpsr=${psr} '{print "/fred/oz002/users/mmiles/MTM_trials/"awkpsr"/for_timing/"$0"*32ch"}')  >> ${psr}.tim_MTM

    cp /fred/oz002/users/mmiles/MPTA_GW/partim_frank/pp_8/${psr}.par .

    touch transfer_finished
fi
