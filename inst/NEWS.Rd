% Check from R:
%  news(db = tools:::.build_news_db_from_package_NEWS_Rd("~/R/Pkgs/Matrix/inst/NEWS.Rd"))
\name{NEWS}
\title{News for \R Package \pkg{asht}}
\encoding{UTF-8}
\newcommand{\CRANpkg}{\href{http://CRAN.R-project.org/package=#1}{\pkg{#1}}}
%% NB: The date (yyyy-mm-dd) is the "Packaged:" date in ../DESCRIPTION


%\section{Version 0.9.7 (???)}{
% \subsection{Changes}{
%  \itemize{
%   \item TO DO: finish poisson2sampleDiff 
%         (see /PUBLISHED/2015 Fay Proschan Brittain/R/poisson2sampleDiff)
%         the meldCD function there is the same that is here.
%  }
% } 
%}



\section{Version 1.0.0 (2023-02-06)}{
 \subsection{Changes}{
  \itemize{
   \item Added WprevSeSp function.
  }
 } 
}

\section{Version 0.9.8 (2021-09-15)}{
 \subsection{Changes}{
  \itemize{
   \item Fixed bug in simulateSS function (see comment inside estNpower function within simulateSS).
  }
 } 
}



\section{Version 0.9.7 (2021-01-25)}{
 \subsection{Changes}{
  \itemize{
   \item Added simulateSS function.
  }
 } 
}



\section{Version 0.9.6 (2020-07-03)}{
 \subsection{Changes}{
  \itemize{
   \item Added unconditional difference option for the statistic for the sign test. This calls ncnemarExactDP, that is why the package now requires exact2x2 (>=1.6.4).
  }
 } 
}




\section{Version 0.9.5.1 (2020-06-22)}{
 \subsection{Changes}{
  \itemize{
   \item added utol argument to meldCD.
   \item Changed wmwTest so that it gives a warning if the uniroot function may fail for method='exact.mc'. See Note section on the man page. 
  }
 } 
}




\section{Version 0.9.4 (2018-09-14)}{
 \subsection{Changes}{
  \itemize{
   \item Added meldtTest and meldCD functions to do two-sample tests by melding.
   \item Added tukeyWelsch function to do multiplicity adjusted pairwise tests from k-samples
   \item Added anovaOneWay, to output an ICC or other estimate, also allow var.eqaul=FALSE using the Brown-Forsythe method.
   \item Added ama1c1cpg data set.
   \item Fixed wmwTest when exact and there are tied values and value is most extreme. 
   For example,  wmwTest(c(1,3),c(3,4))  previously had given an error for calculating 
   the confidence interval. 
   \item use match.arg with method in wmwTest, so that method='asy' will work for method='asymptotic'
   \item updated wmwTest reference
   \item added output to bfTest so it outputs the R statistic.
  }
 } 
}


\section{Version 0.9.3 (2017-10-06)}{
 \subsection{Changes}{
  \itemize{
   \item Added wmwTest to calculate the Wilcoxon-Mann-Whitney test with confidence intervals for the Mann-Whitney parameter that are unified with the test.
  }
 } 
}


\section{Version 0.9.2 (2017-08-10)}{
 \subsection{Changes}{
  \itemize{
   \item Added the abcnonHtest function to give nonparametric ABC (approximate bootstrap confidence) intervals with p-values and estimates.
  }
 } 
}



\section{Version 0.9.1 (2017-05-01)}{
 \subsection{Changes}{
  \itemize{
   \item Fixed error in one-sided tests for wspoissonTest when midp=TRUE. Thanks to Michael Nelson for pointing out this error. Also fixed data description for wspoissonTest.
  }
 } 
}



\section{Version 0.9 (2017-04-01)}{
 \subsection{Changes}{
  \itemize{
   \item Created bfTest function to perform the Behrens-Fisher test using the Behrens-Fisher solution with p-values and CIs.  
  }
 } 
}



\section{Version 0.8 (2017-02-06)}{
 \subsection{Changes}{
  \itemize{
   \item Made changes to the estimate and confidence interval for  \code{wsrTest}. 
   The default estimate previously used the estimator of the pseudo-median from the difference using 
   the estimator from \code{wilcox.test}. 
   The new estimator is based on the shift value where the one-sided p-values are equal. This ensures that the 
   estimate will be within the confidence interval and will not change as the null hypothesis changes.
   \item Updated citation for Fay and Kim (2017).
  }
 } 
}




\section{Version 0.7 (2016-11-15)}{
 \subsection{Changes}{
  \itemize{
   \item Corrected error in \code{wsrTest}. The p-value was incorrectly calculated when mu was not 0, now it is correct. Also added the tieDigits option.
   \item Updated citation for Fay and Kim (2017).
  }
 } 
}





\section{Version 0.6 (2016-05-04)}{
 \subsection{Changes}{
  \itemize{
   \item Modified \code{wspoissonTest} to be calculated using uniroot by default, instead of by Monte Carlo simulation.
  }
 } 
}


\section{Version 0.5 (2015-06-26)}{
   First CRAN release of \pkg{asht}, \dQuote{Applied Statistical Hypothesis Tests} by Michael P. Fay (author, maintainer).
  \subsection{Tests}{
    \itemize{
      \item    \code{signTest} gives exact sign test with confidence intervals
      \item \code{cvTest} gives coefficient of variance tests based on normality or lognormal assumptions
     \item \code{quantileTest} and \code{medianTest} give exact tests on a quantile (including median). No assumptions on the distribution are needed (ordered factors are allowed). 
     \item \code{var1Test} gives one-sample test of normal variance. 
    \item \code{wsrTest} gives exact Wilcoxon signed rank test, gives exact matching confidence intervals.
    \item \code{wspoissonTest} gives confidence intervals for weighted sums of Poisson, useful for directly standardized rates.      
    }
  }
}
