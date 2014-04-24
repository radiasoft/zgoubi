      character txt80*80, txtfig*20, txtM*4, txtQ*20
      parameter (lunW=8)
      integer debstr, finstr

      open(unit=lunW,file='./Log/log.tex')

      txt80='\documentclass[11pt]{article}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\usepackage{draftcopy}'
      write(lunW,*) txt80
      txt80='%\usepackage[draft]{graphicx}'
      write(lunW,*) txt80
      txt80='\usepackage{graphicx}'
      write(lunW,*) txt80
      txt80='\usepackage{wrapfig}'
      write(lunW,*) txt80
      txt80='\usepackage{amssymb}'
      write(lunW,*) txt80
      txt80='\usepackage{lscape}'
      write(lunW,*) txt80
      txt80='\usepackage{times}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\oddsidemargin =-0.1in'
      write(lunW,*) txt80
      txt80='\evensidemargin=-0.1in'
      write(lunW,*) txt80
      txt80='\textwidth=6.8in'
      write(lunW,*) txt80
      txt80='\textheight=10.2in'
      write(lunW,*) txt80
      txt80='\topmargin=-.8in'
      write(lunW,*) txt80
      txt80='\footskip=0in'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\newcommand{\Br}{\ensuremath{B\!\rho}}'
      write(lunW,*) txt80
      txt80='\newcommand{\bull}{\ensuremath{\bullet~}}'
      write(lunW,*) txt80
      txt80='\newcommand{\nib}{\noindent \ensuremath{\bullet~}}'
      write(lunW,*) txt80
      txt80='\newcommand{\com}{{center of mass}}'
      write(lunW,*) txt80
      txt80='\newcommand{\hbrk}{\hfill \break}'
      write(lunW,*) txt80
      txt80='\newcommand{\lab}{{laboratory frame}}'
      write(lunW,*) txt80
      txt80='\newcommand{\MC}{Monte~Carlo}'
      write(lunW,*) txt80
      txt80='\newcommand{\rms}{\ensuremath{rms}}'
      write(lunW,*) txt80
      txt80='\newcommand{\wrt}{{with respect to}}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\newcommand{\rfrncA}{\rm  }'
      write(lunW,*) txt80
      txt80='\newcommand{\rfrncB}{\rm }'
      write(lunW,*) txt80
      txt80='\newcommand{\rfrncC}{\rm }'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\pagestyle{myheadings}'
      write(lunW,*) txt80
      txt80='\markboth{\small \rfrncA ~ ~ \rfrncB ~ ~  \rfrncC \hfill}'
      write(lunW,*) txt80
      txt80='         {\small \rfrncA ~ ~ \rfrncB ~ ~  \rfrncC \hfill}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\pagestyle{headings}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\begin{document}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\thispagestyle{myheadings}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='~ ~ ~ '
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80=' ~ ~ ~ \hfill F. M\''eot, BNL, Sept. 2009'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\vfill'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\begin{center}'
      write(lunW,*) txt80
      txt80='\bf \large '
      write(lunW,*) txt80
      txt80='Spin tracking in AGS based on ray-tracing methods'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='- bare lattice, no snakes -'
      write(lunW,*) txt80
      txt80='\end{center}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\vfill'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\tableofcontents'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\vfill'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\clearpage'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\section{Introduction}'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='      end'
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80=''
      write(lunW,*) txt80
      txt80='\clearpage'
      
      close(lunW)
      stop
      end
