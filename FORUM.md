--------------------------------------------------------------------------------

Welcome to the match!

Here is the competitor pack

Inside the pack input and ground-truth folders contain training datasets and the corresponding ground truths. Inside scorer folder there are sources for Docker scorer container, and run.sh bash script runs the testing and scoring locally, using cmap/sig_dpeak_tool Docker container as a benchmark solution (notice that after pulling that container is tagged as cmap/solution, and it is reference under that name further in the script).

Notice, that we expect to change the way the final score is aggregated from COR, AUC, and peformance components. We will announce the exact formula, and provide updated scoring script in the next few days.

At the first run run.sh will spend some time pulling the benchmark solution container, and building the scorer from sources; on subsequent runs it will just reuse them. The benchmark solution is based on a stripped down version of DPeak algorithm, which runs significantly faster than 30 minutes on 12-cores, mentioned in the challenge description, but produces somewhat less accurate results.

To submit into this challenge, you will need to dockerize your solution, to expose the same interface as cmap/sig_dpeak_tool container exposes, and you will submit ZIPped source code of your Docker image to Topcoder system. If you are not familiar with Docker, have a look at the sources of scorer container and Docker documentation. If you have any difficulties to get it work, I'll be happy to help, advicing on the dockerization part.

Here I need to apologies that we don't have the online scoring setup and leaderboard live at the moment of starting the match. It turns out that we need some extra efforts to deploy the scorer to Topcoder system, and we plan to get it ready early next week (in case of delay, we will extend the match accordingly). However, we are sure about the problem itself, thus we decided to launch the match already. Once the online scoring system is up, we will announce it in the forum, along with the prize eligibility threshold, which will be a bit above the benchmark score (presumably ~10% higher, which should be relatively easy to achieve).

The last, and not the least, notice that multithreading is allowed in this match. We plan to run provisional scoring on 8-cores machines. The exact specifications will be announced once the online scoring is up.

Good luck! :)

--------------------------------------------------------------------------------

There won't be a training phase. If you try an algorithm relying on machine learning, the model should be trained locally, and all necessary model assets should be packed into the dockerized submission, so that they can be just loaded at execution time. Then, if you win a prize, the (re-)training procedure should be documented in the manual submitted after the match as a part of the final deliverable.

--------------------------------------------------------------------------------

After some additional discussions between the TC and client teams, along with some testing and fixes on the way, I am happy to present the final formula for the overal score aggregation:

SCORE = MAX_SCORE * (MAX(COR, 0))^2 * AUC^2 * exp(- T_run / (3 * T_benchmark) )

where:

MAX(COR, 0) is just an inline way to write: COR, if COR positive; zero, if COR is negative;
T_benchmark is the average test case deconvolution time spent by the benchmark solution. Two benchmark time values (one per test case) are hardcoded at the lines ##18,19 of scorer/aggregate_score.R script in the updated competitor pack. For the purposes of local scoring on your machine, you can replace thes bencmark times by the actual values you get on your machine. When you run the testing & scoring code, the times are written out to the last columns in the output/scoring_dpk.csv and output/scoring_litmus.csv files.


The individual values of AUC and COR calculated by the scoring script also changed a bit, as we fixed some bugs in the initial script version.


Here I was planning to write about the prize elibility threshold score. Hovewer, I found, that benchmark scores with training, provisional, and final testing data are noticeably different, thus proclaiming a single threshold score makes no sense. Let's put it this way: to be eligible for a prize, your algorithm should demostrate a clear improvement over the benchmark, so that just reproducing the benchmark solution won't make you eligible. In practice, we expect that it will be relatively easy to beat the benchmark speed (as the benchmark code is implemented in MatLab), thus many competitors will show significant score improvements, and we won't have to bother about the exact elibility threshold. In case we have to figure it out, an algorithm 10% faster than bencmark, and with the same accuracy, will be considered eligible.


Please, find the updated competitor pack here: https://drive.google.com/open?id=1865JSqF2NqxOcEMcDV57a10k6YTWRm6m
If you want to save some traffic, here is the light version (does not contain inputs and ground truth files, which were not changed): https://drive.google.com/open?id=16Xf4McmIfxMyMR_dzVqEBLdR5kZ7o2hs


Once again, my apologies for not providing this at the challenge launch time.

The online scoring system and the leaderboard should be, hopefully, ready soon. Stay tuned for further announcements.

