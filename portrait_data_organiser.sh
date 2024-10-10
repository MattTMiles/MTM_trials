#!/bin/bash

psr=$1

#Collect data

if [[ ! -f "${psr}rm_trials/transfer_finished" ]]; then

    # echo $psr 
    # for dir in $(ls -d /fred/oz005/timing_processed/$psr/20*); do 
    #     if [[ -f $(echo $dir/*/decimated/*zap*928ch_4p_1t.ar) ]]; then 
    #         cp $dir/*/decimated/*zap*928ch_4p_1t.ar $psr/rm_trials/data/; 
    #     elif [[ -f $(echo $dir/*/decimated/*zap.*928ch_4p_1t.ar) ]]; then
    #         cp $dir/*/decimated/*zap*928ch_4p_1t.ar $psr/rm_trials/data/;
    #     else
    #         echo "${dir} doesn't contain file"
    #         echo "${dir}" >> ${psr}/${psr}_missing.list
    #     fi; 
    # done; 
    # #done


    cd ${psr}/rm_trials/data

    #Change data into correct format and dlyfix

    # for arch in $(ls *ar);
    # do
    #     echo ${arch}
    #     pam -mpT ${arch}
    #     /fred/oz002/users/mmiles/dlyfix/dlyfix ${arch} -u . -o ${arch%.ar}.dly
    # done

    #rm J*ar
    #Apply the correct dispersion measure to the data, also apply the RM
    for arch in $(ls J*ar); 
    do 
        echo ${arch} 
        pdmp -g /ps -bf ${arch} 
        dm_temp=$(cat pdmp.best | head -n 9 |tail -n 1 | awk '{print($1)}')
        pam -md ${dm_temp} ${arch} --update_dm
        pam -D ${arch} -e dly_D
    done

    for arch in $(ls J*dly_D); 
    do
        rmfit -L -r ${arch} #output as rmfit.out
        rm_temp=$(cat rmfit.out | awk '{print $12}' | tail -n1)
        pam --aux_rm ${rm_temp} -e ar_aux_rm ${arch}
        
        pam --RM $rm_temp -e ar_rm ${arch}
        pam -mR $rm_temp ${arch%.dly_D}.ar_rm

        rm rm*out #get rid of files before next iteration
        rm *ps

    done

    #Need to change the keyword for the other RM attempts
    seed_dir="/fred/oz002/users/mmiles/MPTA_GW/partim_august23_snr10/partim_updated/"
    psradd -PT $(psrstat -Q -c snr $(awk '{print $1}' ${seed_dir}/${psr}.tim | uniq | cut -d "/" -f 9 | grep J | rev | cut -f 2- -d '.' | rev | sed 's/......$//' | awk -v awkpsr=${psr} '{print "*"$0"*.ar_aux_rm"}') | sort -g -k 2 | tail -n $(echo $(( $(ls J*.ar_aux_rm| wc -l) /1))) | awk '{print $1}') -o ${psr}_grand_T.ar_aux_rm



    #Create grand data object

    #seed_dir="/fred/oz002/users/mmiles/MPTA_GW/partim_august23_snr10/partim_updated/"
    #psradd -PT $(psrstat -Q -c snr $(awk '{print $1}' ${seed_dir}/${psr}.tim | uniq | cut -d "/" -f 9 | grep J | rev | cut -f 2- -d '.' | rev | awk -v awkpsr=${psr} '{print "*"$0"*.dly_D"}') | sort -g -k 2 | tail -n $(echo $(( $(ls J*fluxcal.dly_D | wc -l) /1))) | awk '{print $1}') -o ${psr}_grand_T.dly_D


fi
