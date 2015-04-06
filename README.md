# camb - chemistry aware model builder
======

camb - chemistry aware model builder is an R package that can be used for the rapid generation of quantitative predictive models in the area of medicinal chemistry (QSAR, QSPR, QSAM, PCM). It is aimed at both advanced and beginner R users.
Its capabilities include the standardisation of representation of chemical structures, computation of 905 two-dimensional and 14 fingerprint type descriptors for small molecules, 8 types of amino acid descriptors, 13 whole protein sequence descriptors, filter methods for feature selection, generation of predictive models (R package {\it caret}), as well as techniques to ensemble these models (R package {\it caretEnsemble}).
Results can be visualised through high-quality, customisable plots (R package {\it ggplot2}).

This is the root folder which holds the package folder as well as other folders which contain examples of package use.

Two tutorials concerning the application of camb in the context of QSPR and Proteochemometrics are available in the examples folder.

Coding is done with the Google's R style guide: http://google-styleguide.googlecode.com/svn/trunk/Rguide.xml#functiondefinition

# INSTALLATION TIPS:

camb can be installed the by typing: library(devtools); install_github("cambDI/camb/camb")

Please read the following operating system specific instructions if the above method fails for you at some point.

1. OSX: You'll need to have the Xcode installed. If you don't already, this install should be triggered automatically during the installation although you may need to restart the installation: install_github("cambDI/camb/camb")
If it does not install automatically, you can install it from here: https://developer.apple.com/xcode/downloads/
relevant error: 'xcode-select: note: no developer tools were found at '/Applications/Xcode.app', requesting install. Choose an option in the dialog to download the command line developer tools.'

You'll also need the boost library installed. This can installed through the brew package manager by typing 'brew install boost' in a terminal window. If you don't have brew install yet you can type 'ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"' to install it. If you already use macports (another package manager) you can install boost by typing 'sudo port install boost' in a terminal window.

2. CentOs: Make sure that cmake (version >= 2.8) is installed and that the package rJava can be installed (in case you want to use the function GeneratePadelDescriptors). CentOs users have reported issues in this regard. The C++ library boost also needs to be installed: http://www.boost.org/
Please email us with any other problems you might encounter! Thanks!


