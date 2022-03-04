FROM rocker/verse:4.1.0

# adding plotting packages to from day 1 demo
RUN install2.r --error \
    --deps TRUE \
    ggrepel \
    ggthemes \
    here

## adding data grabbing packages
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    wget curl git 

RUN install2.r --error --deps TRUE \
    tableone \
    ggridges 


RUN install2.r --deps TRUE \    
    tidygraph 

RUN install2.r --error --deps TRUE \
    ggraph 

RUN install2.r --error --deps TRUE \        
    cowplot  

## install stuff at the end of Dan's script
RUN install2.r --error --deps TRUE \
    tidymodels \
    vip \
    kernlab

## install the other odd things
RUN install2.r --error --deps TRUE \
    MatchIt \
    effsize