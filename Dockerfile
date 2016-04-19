
########################################################
# Dockerfile for automated build of latest SONAR code  # 
# Uses the scharch/sonar-base Docker image, which in   # 
#  turn is based on Ubuntu with various SONAR          #
#  dependecies installed.                              #
# Please see https://github.com/scharch/SONAR for more #
#  information.                                        #
########################################################

FROM scharch/sonar-base:latest
MAINTAINER Chaim Schramm cs3037@columbia.edu

#pull latest SONAR source code
WORKDIR /sonar
RUN git pull https://github.com/scharch/SONAR.git
