# Detecting-changes-in-ConditionalCorrelations
This simple package provides the Matlab code to implement the tests for detecting change points in conditional correlation structure.
Reference:
Barassi, M., Horváth, L. and Zhao, Y., (2018). Change‐Point Detection in the Conditional Correlation Structure of Multivariate Volatility Models. Journal of Business and Economic Statistics. **38**, 340-349.
The package contains two main files, cv_stat.m and uyility.m

## cv_stat.m
The file provides codes for computing the test statistics and the critical values for assessing the null hypothesis of no change.

## utility.m
The file provides all codes that need to run the main functions in the cv_stat.m file. Note that the MFE toolbox is also required to run this code, and it is publically accessible from the link "https://www.kevinsheppard.com/MFE_Toolbox".
