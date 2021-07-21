## Description
This repository includes all the code I have done in *2017* for my Bachelor's Degree Final project. The project consisted of my own implementation in MATLAB of the algorithm for unsupervised source separation based on the concept of Average Harmonic Structure (AHS) modelling that has been proposed on the following paper:

* Z. Duan, Y. Zhang, C. Zhang and Z. Shi, "__Unsupervised Single-Channel Music Source Separation by Average Harmonic Structure Modeling__", in _IEEE Transactions on Audio, Speech, and Language Processing_, vol. 16, no. 4, pp. 766-778, May 2008, [doi: 10.1109/TASL.2008.919073](https://doi.org/10.1109/TASL.2008.919073).

The results I obtained using my own implementation were quite similar to the results of the authors. I have written a complete report (in Portuguese) of this project, which can be found [here](http://monografias.poli.ufrj.br/monografias/monopoli10022740.pdf)

## Previous Code Utilised
Since the original authors did not provide the original code they utilised in their experiments, I have done my own implementation and used the same audio signals that they used in the paper in my tests. The original audio files can be downloaded from the [website with the paper results](http://mperesult.googlepages.com/musicseparationresults).

Moreover, part of the code (mostly the part related to AHS computation) has been taken from another publication, which also used the concept of AHS, but focused on multi-pitch estimation. 

* Z. Duan, B. Pardo, and C. Zhang, "__Multiple Fundamental Frequency Estimation by Modeling Spectral Peaks and Non-Peak Regions__", in _IEEE Transactions on Audio, Speech, and Language Processing_, vol. 18, no. 8, pp. 2121-2133, 2010

The full code related to this publication can be found on [Zhiyao Duan's personal webpage](http://www2.ece.rochester.edu/projects/air/publications.html). I have taken some of their code, adapted for this project and added in this repo. In the header of each file or function you can find the information regarding its authorship if it has not been originally created by me. 

## Main Usage
The main code is the file `MainSeparationAHS.m`, which will read the parameters of the configuration file `SetParameters_music.m` and will perform the whole analysis of the mixture file. The mixture file should be set up by editing the macro `AUDIO_FILE` in `line 14` of the script `SetParameters_music.m`. As default, this line reads

```
AUDIO_FILE = 'Mixture.wav'
```

Detailed explanation of the other parameters in this script can be found in the code as comments.