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
    printf "\n\t%-9s  %-40s"  "0.2.1"     ""
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
           ./bin/exe pwa
	       ;;

esac

}

sub_0_2() {

case $option in
    
    # --------------------------------------------------------------------------
    #   
    # --------------------------------------------------------------------------

    0.2.1) echo "..."
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
    #  
    # --------------------------------------------------------------------------

    0.2) echo "..."
         usage_0_2 
         echo "Please enter your option: " 
         read option  
         sub_0_2 option 
	     ;;

    0.2.*) echo "..."
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
