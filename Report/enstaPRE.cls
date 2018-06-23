\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{enstaPRE}
  [2014/05/12 - v1 - Rapport de PRE - Quentin Byron - ENSTA ParisTech, Promotion 2015]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        CLASS  LOADER        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\LoadClass[a4paper,12pt,twoside,openright]{report}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%      REQUIRED PACKAGES      %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\RequirePackage{array}
\RequirePackage[english]{babel}
\RequirePackage{datatool}
\RequirePackage{environ}
\RequirePackage[T1]{fontenc}
\RequirePackage{fancyhdr}
\RequirePackage[top=3cm, bottom=3cm, left=3cm, right=3cm]{geometry}
\RequirePackage{graphicx}
\RequirePackage[
    bookmarks         = true,%     % Signets
    bookmarksnumbered = true,%     % Signets numérotés%
    bookmarksopen     = true,%     % Signets ouverts
    colorlinks        = true,%     % Liens en couleur : true ou false
    urlcolor          = blue,%     % Couleur des liens externes
    linkcolor         = black,%    % Couleur des liens internes
    citecolor         = black,%    % Couleur des citations
    ]{hyperref}% 
\RequirePackage{ifthen}
\RequirePackage{indentfirst}
\RequirePackage[toc,nopostdot,nonumberlist]{glossaries}
\RequirePackage[utf8]{inputenc}
\RequirePackage{stmaryrd}

\def\enstaLang{en}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        LANGUAGE  DEF        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\DeclareOption{francais}{\def\enstaLang{fr} \PassOptionsToPackage{francais}{babel}}
\DeclareOption{english}{\def\enstaLang{en} \PassOptionsToPackage{english}{babel}}
\ProcessOptions

\ifthenelse{\equal{\enstaLang}{fr}}{
  \addto\captionsenglish{\renewcommand{\chaptername}{Partie}}
  \def\enstaLangPRE{Projet de Recherche (PRE)}
  \def\enstaLangSpecialite{Spécialité :}
  \def\enstaLangAnneeScolaire{Année Scolaire :}
  \def\enstaLangMentionConfidentialite{Mention de confidentialité}
  \def\enstaLangAuteur{Auteur :}
  \def\enstaLangTuteurENSTA{Tuteur ENSTA ParisTech :}
  \def\enstaLangPromo{Promotion :}
  \def\enstaLangTuteurOrganisme{Tuteur organisme d'accueil :}
  \def\enstaLangStageEffectueDu{Stage effectué du}
  \def\enstaLangStageEffectueAu{au}
  \def\enstaLangNomOrganisme{Nom de l'organisme d'accueil :}
  \def\enstaLangAdresse{Adresse :}
}{
  \addto\captionsenglish{\renewcommand{\chaptername}{Part}}
  \def\enstaLangPRE{Research Internship (PRE)}
  \def\enstaLangSpecialite{Field of Study:}
  \def\enstaLangAnneeScolaire{Scholar Year:}
  \def\enstaLangMentionConfidentialite{Confidentiality Notice}
  \def\enstaLangAuteur{Author:}
  \def\enstaLangTuteurENSTA{ENSTA ParisTech Tutor:}
  \def\enstaLangPromo{Promotion:}
  \def\enstaLangTuteurOrganisme{Host Organism Tutors:}
  \def\enstaLangStageEffectueDu{Internship from}
  \def\enstaLangStageEffectueAu{to}
  \def\enstaLangNomOrganisme{Name of the host organism:}
  \def\enstaLangAdresse{Address:}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        DEFAULT  FONT        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand*\familydefault{\rmdefault}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%          VARIABLES          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\logo}[1]{\def\enstaLogoCorp{#1}}
\newcommand{\specialite}[1]{\def\enstaSpecialite{#1}}
\newcommand{\annees}[1]{\def\enstaAnneeScol{#1}}
\newcommand{\titre}[1]{\def\enstaTitre{#1}}
\newcommand{\soustitre}[1]{\def\enstaSousTitre{#1}}
\newcommand{\illustration}[1]{\def\enstaIllustration{#1}}
\newcommand{\confidentialite}[1]{\def\enstaConfidentialite{#1}}
\newcommand{\auteur}[1]{\def\enstaAuteur{#1}}
\newcommand{\promotion}[1]{\def\enstaPromotion{#1}}
\newcommand{\tuteurENSTA}[1]{\def\enstaTuteurENSTA{#1}}
\newcommand{\tuteurOrganisme}[1]{\def\enstaTuteurOrganisme{#1}}
\newcommand{\dateDebut}[1]{\def\enstaDateDebut{#1}}
\newcommand{\dateFin}[1]{\def\enstaDateFin{#1}}
\newcommand{\organisme}[1]{\def\enstaOrganisme{#1}}
\newcommand{\adresseOrganisme}[1]{\def\enstaAdresseOrganisme{#1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        HEADER/FOOTER        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\pagestyle{fancy}
\fancyhead{}
\fancyfoot{}
\fancyhead[CO]{\leftmark}
\fancyhead[CE]{\enstaTitre}
\fancyfoot[C]{\enstaAuteur ~/ \enstaOrganisme \\ \textbf{\color{red}\enstaConfidentialite}}
\fancyfoot[RO, LE] {\thepage}
\fancypagestyle{plain}{
\fancyhead{}
\fancyfoot{}
\fancyhead[CO]{\enstaTitre}
\fancyfoot[C]{\enstaAuteur ~/ \enstaOrganisme \\ \textbf{\color{red}\enstaConfidentialite}}
\fancyfoot[RO, LE] {\thepage}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           CHAPTER           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\renewcommand\thechapter{\Roman{chapter}}
\newcommand{\partb}[1]{\chapter*{#1\addcontentsline{toc}{chapter}{#1}\markboth{#1}{#1}}}
\renewcommand{\part}{\chapter}

\def\cleardoublepage{\clearpage
 \if@twoside
  \ifodd\c@page\else
   \null\thispagestyle{chapterverso}\newpage
   \if@twocolumn\null\newpage\fi
   \fi
  \fi
 }%
\def\ps@chapterverso{\ps@empty}%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         FIRST  PAGE         %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\couverture}{
  \begin{titlepage}
  {\bfseries 
    \begin{center}
    \vspace*{\stretch{2}}
    \begin{tabular*}{\textwidth}{@{\extracolsep{\fill}} lr}
      \includegraphics[height=4cm]{logos/logo-ensta.png} & \includegraphics[height=4cm]{\enstaLogoCorp}
    \end{tabular*} \\
    \vspace*{\stretch{0.5}}
    \Huge{\enstaLangPRE} \\
    \vspace*{\stretch{0.5}}
    \Large \enstaLangSpecialite ~\enstaSpecialite \\
    \enstaLangAnneeScolaire ~\enstaAnneeScol \\
    \vspace*{\stretch{0.8}}
    \Huge \enstaTitre \\
    \vspace*{\stretch{0.5}}
    \includegraphics[height=4cm]{\enstaIllustration} \\
    \vspace*{\stretch{0.5}}
    {\color{red}\enstaLangMentionConfidentialite \\
    \Large \enstaConfidentialite} \\
    \vspace*{\stretch{0.5}}
    \begin{minipage}{0.49\textwidth}
      \begin{center}
  \Large \enstaLangAuteur \\
  \enstaAuteur\\
  \vspace*{0.3cm}
  \large \enstaLangTuteurENSTA\\
  \enstaTuteurENSTA \\
      \end{center}
    \end{minipage}
    \begin{minipage}{0.49\textwidth}
      \begin{center}
  \Large \enstaLangPromo \\
  \enstaPromotion \\
  \vspace*{0.3cm}
  \large \enstaLangTuteurOrganisme \\
  \enstaTuteurOrganisme \\
      \end{center}
    \end{minipage} \\
    \vspace*{\stretch{1}}
    \Large \enstaLangStageEffectueDu ~\enstaDateDebut ~\enstaLangStageEffectueAu ~\enstaDateFin \\
    \vspace*{\stretch{0.5}}
    \large
    \begin{tabular}{l}
      \begin{tabular}{l@{}l}
       \enstaLangNomOrganisme & \begin{tabular}{l}\enstaOrganisme\end{tabular}
      \end{tabular} \\
      \begin{tabular}{l@{}l}
       \enstaLangAdresse & \begin{tabular}[t]{l}\enstaAdresseOrganisme\end{tabular}
      \end{tabular}
    \end{tabular}
    \vspace*{\stretch{2}}
    \end{center} }
  \end{titlepage}
  \newgeometry{top=2.5cm, bottom=2.5cm, left=2.5cm, right=2.5cm}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   CONFIDENTIALITY  NOTICE   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\pageConfidentialite}[1]{\partb{\enstaLangMentionConfidentialite}\vfill#1\vfill}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           GLOSSARY          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\afficheGlossaire}{\printglossaries\cleardoublepage}
\newcommand{\glossaire}[1]{\makeglossaries\include{#1}}
\newcommand{\ac}[1]{\acrshort{#1}\index{#1@\gls{#1}}}
\newcommand{\gl}[1]{\gls{#1}\index{#1@\gls{#1}}}
\newcommand{\Gl}[1]{\Gls{#1}\index{#1@\gls{#1}}}
\newcommand{\glp}[1]{\glspl{#1}\index{#1@\gls{#1}}}
\newcommand{\Glp}[1]{\Glspl{#1}\index{#1@\gls{#1}}}
\newcommand{\gll}[2]{\gllink{#1}{#2}\index{#1@\gls{#1}}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%         BIBLIOGRAPHY        %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\bibliographie}[1]{\bibliographystyle{unsrt}\bibliography{#1}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            INDEX            %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\makeindex
\newcommand{\afficheIndex}{\printindex}
\newcommand{\id}[1]{\index{#1}#1}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    LIST OF TABLES/FIGURES   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcommand{\listedestables}{\listoftables}
\newcommand{\listedesfigures}{\listoffigures}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           TABULAR           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\newcolumntype{C}[1]{>{\centering\arraybackslash}m{#1}}
% \begin{tableau}{<label>}{<name>}{<columns>}
\NewEnviron{tableau}[3]
{\begin{table}[!h]\centering\begin{tabular}{#3}
  \BODY
\end{tabular}\caption{\label{tab:#1} #2}\end{table}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            FIGURE           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%\image{<label>}{<name>}{<options>}{<filename>}
%\imageDeux{<label1>}{<name1>}{<options1>}{<filename1>}{<label2>}{<name2>}{<options2>}{<filename2>}
\newcommand{\image}[4]{\begin{figure}[!h]\centering\includegraphics[#3]{#4}\caption{\label{fig:#1}#2}\end{figure}}
\newcommand{\imageDeux}[8]{\begin{figure}[!h]\begin{minipage}{0.49\linewidth}\centering\includegraphics[#3]{#4}\caption{\label{fig:#1}#2}\end{minipage}\hfill\begin{minipage}{0.49\linewidth}\centering\includegraphics[#7]{#8}\caption{\label{fig:#5}#6}\end{minipage}\end{figure}}
\NewEnviron{dessin}[2]{
\begin{figure}[!h]\centering\begin{pspicture}
  \BODY
\end{pspicture}\caption{\label{fig:#1} #2}\end{figure}
}
