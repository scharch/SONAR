
########################################################
# Dockerfile for automated build of latest SONAR code  # 
# Uses the scharch/sonar-base Docker image, which in   # 
#  turn is based on Ubuntu with various SONAR          #
#  dependecies installed.                              #
# Please see https://github.com/scharch/SONAR for more #
#  information.                                        #
########################################################

FROM biopython/biopython
MAINTAINER Chaim Schramm chaim.schramm@nih.gov

#add docopt
RUN pip3 install docopt

#install libraries for bioperl
RUN apt-get update && apt-get install -y \
  build-essential \
  gcc-multilib \
  perl \
  cpanminus \
  libdb-dev \
  graphviz 

#install perl modules that are prerequisites
RUN cpanm \
  CPAN::Meta \
  YAML \
  Digest::SHA \
  Module::Build \
  Test::Most \
  Test::Weaken \
  Test::Memory::Cycle \
  Clone \
  HTML::TableExtract \
  Algorithm::Munkres \
  Algorithm::Combinatorics \
  Statistics::Basic \
  List::Util \
  Array::Compare \
  Convert::Binary::C \
  Error \
  Graph \
  GraphViz \
  Inline::C \
  PostScript::TextBlock \
  Set::Scalar \
  Sort::Naturally \
  Math::Random \
  Spreadsheet::ParseExcel \
  IO::String \
  JSON \
  Data::Stag \
  CGI \
  Bio::Phylo

#now actually install BioPerl
RUN cpanm -v \
    https://github.com/bioperl/bioperl-live/archive/master.tar.gz

#install PyQt and ete3
RUN apt-get install -y \
    python3-pyqt4 python3-pyqt4.qtopengl python-lxml python-six

RUN pip3 install --upgrade ete3

#install R
RUN echo "deb http://cran.rstudio.com/bin/linux/ubuntu xenial/" >> /etc/apt/sources.list
RUN gpg --keyserver keyserver.ubuntu.com --recv-key E084DAB9
RUN gpg -a --export E084DAB9 | apt-key add -
RUN apt-get update && apt-get install -y r-base r-base-dev

#install R packages
RUN R --vanilla -e 'install.packages(c("docopt","MASS","ggplot2"), repos="http://cran.us.r-project.org")'


#pull latest SONAR source code
WORKDIR /sonar
RUN git clone https://github.com/scharch/SONAR.git
