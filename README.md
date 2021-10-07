# Parallel Smith-Waterman with OpenMP
This is a parallel Smith-Waterman algorithm written for COMP90025. 

Interesting notes:
  * NUMA aware memory management ensures NUMA costs are mitigated
  * 13x speedup over sequential algorithm


Grade:
  - Code: 10/10
  - [Report](https://github.com/isubasinghe/parallel-smith-waterman/blob/main/code/isithasubasinghe_report.pdf) 9.25/10

Ranking:
  No official rankings were published but I believe this is the **fastest implementation in the cohort**. 
  My speed test put the timing of this at 108 seconds and the fastest was reported to be "around 105 seconds", this deviation seems enough that it is plausible that this is my implementation. 
  
  This is also the highest grade in the class for COMP90025 at 19.25/20. This was **verified** by checking the histogram of results, only one made it into the 19-20 bucket which leads
  to the natural conclusion that I must have received the highest grade for this assignment.
