#!/usr/bin/env bash

# Main driver to submit jobs 
# Author SHI Xin <shixin@ihep.ac.cn>
# Modified by JING Maoqiang <jingmq@ihep.ac.cn>
# Created [2020-05-29 May 09:58]

usage() {
    printf "NAME\n\tsubmit.sh - Main driver to submit jobs\n"
    printf "\nSYNOPSIS\n"
    printf "\n\t%-5s\n" "./submit.sh [OPTION]" 
    printf "\nOPTIONS\n" 
    printf "\n\t%-9s  %-40s"  "0.1"       "[setup PWA]"
    printf "\n\t%-9s  %-40s"  "0.2"       "[execute PWA]"
    printf "\n\n" 
}

usage_0_1() {
    printf "\n\t%-9s  %-40s"  ""          ""
    printf "\n\t%-9s  %-40s"  "0.1.1"     "Build pwa executor"
    printf "\n\t%-9s  %-40s"  "0.1.2"     "Prepare four momentum files for data, PHSP MC and background"
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

usage_0_2() {
    printf "\n\t%-9s  %-40s"  ""          ""
    printf "\n\t%-9s  %-40s"  "0.2.1"     "[STEP 1] Calculate fun(i,j) using MC sample"
    printf "\n\t%-9s  %-40s"  "0.2.2"     "[STEP 2] Optimize fractions of resonances"
    printf "\n\t%-9s  %-40s"  "0.2.3"     "[STEP 3] Output weight of MC sample"
    printf "\n\t%-9s  %-40s"  "0.2.4"     "[STEP 4] Tablize the numberical results"
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

usage_0_3() {
    printf "\n\t%-9s  %-40s"  ""          ""
    printf "\n\t%-9s  %-40s"  "0.3.1"     ""
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

usage_0_4() {
    printf "\n\t%-9s  %-40s"  ""          ""
    printf "\n\t%-9s  %-40s"  "0.4.1"     ""
    printf "\n\t%-9s  %-40s"  ""           ""
    printf "\n"
}

if [[ $# -eq 0 ]]; then
    usage
    echo "Please enter your option: "
    read option
else
    option=$1
fi

sub_0_1() {

case $option in
    
    # --------------------------------------------------------------------------
    #  setup PWA 
    # --------------------------------------------------------------------------

    0.1.1) echo "Build pwa executor..."
           ./bin/exe pwa
	       ;;

    0.1.2) echo "Prepare four momentum files for data, PHSP MC and background..."
           cd python
           ./get_momentum.py 4680 Dplus D data
           ./get_momentum.py 4680 Dplus D background
           ./get_momentum.py 4680 Dplus D mc
	       ;;

esac

}

sub_0_2() {

case $option in
    
    # --------------------------------------------------------------------------
    #  execute PWA 
    # --------------------------------------------------------------------------

    0.2.1) echo "[STEP 1] Calculate fun(i,j) using MC sample..."
           cd python
           ./execute_pwa.py step1
	       ;;

    0.2.2) echo "[STEP 2] Optimize fractions of resonances..."
           cd python
           ./execute_pwa.py step2
	       ;;

    0.2.3) echo "[STEP 3] Output weight of MC sample..."
           cd python
           ./execute_pwa.py step3
	       ;;

    0.2.4) echo "[STEP 4] Tablize the numberical results..."
           cd python
           ./execute_pwa.py step4
	       ;;

esac

}

sub_0_3() {

case $option in
    
    # --------------------------------------------------------------------------
    # 
    # --------------------------------------------------------------------------

    0.3.1) echo "..."
	       ;;

esac

}

sub_0_4() {

case $option in
    
    # --------------------------------------------------------------------------
    #   
    # --------------------------------------------------------------------------

    0.4.1) echo "..."
	       ;;

esac

}


case $option in
   
    # --------------------------------------------------------------------------
    #  Setup PWA 
    # --------------------------------------------------------------------------

    0.1) echo "Setting up PWA..."
         usage_0_1 
         echo "Please enter your option: " 
         read option  
         sub_0_1 option 
	     ;;

    0.1.*) echo "Setting up PWA..."
           sub_0_1 option  
           ;;  
        
    # --------------------------------------------------------------------------
    #  Execute PWA
    # --------------------------------------------------------------------------

    0.2) echo "Executing PWA..."
         usage_0_2 
         echo "Please enter your option: " 
         read option  
         sub_0_2 option 
	     ;;

    0.2.*) echo "Executing PWA..."
           sub_0_2 option  
           ;;  

    # --------------------------------------------------------------------------
    #  
    # --------------------------------------------------------------------------

    0.3) echo "..."
         usage_0_3 
         echo "Please enter your option: " 
         read option  
         sub_0_3 option 
	     ;;

    0.3.*) echo "..."
           sub_0_3 option  
           ;;  

    # --------------------------------------------------------------------------
    #  
    # --------------------------------------------------------------------------

    0.4) echo "..."
         usage_0_4 
         echo "Please enter your option: " 
         read option  
         sub_0_4 option 
	     ;;

    0.4.*) echo "..."
           sub_0_4 option  
           ;;  

esac
