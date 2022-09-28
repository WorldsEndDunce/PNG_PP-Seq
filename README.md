# PP-Almost-Seq
This forked repository uses a modified PP-Seq algorithm called "PP-Almost-Seq" to detect repeated firings of polychronous neuronal groups (PNGs) in data and detect the next-best-candidate sequences. PNGs are collections of neurons that fire with a specific order and delays and are theorized to store memories in the brain. Currently, the code focuses on stress-testing the algorithm on multiple instances of planted sequences. Research endeavors are documented [here](https://docs.google.com/document/d/135YVOKIUejlhzhcud6MMVL1u-tEZETtClIymirDUmZI/edit?usp=sharing). To learn more about my research, you can also peruse the slides at the top of the documentation or check out this poster:
<div style="position: relative; width: 100%; height: 0; padding-top: 66.6667%;
 padding-bottom: 0; box-shadow: 0 2px 8px 0 rgba(63,69,81,0.16); margin-top: 1.6em; margin-bottom: 0.9em; overflow: hidden;
 border-radius: 8px; will-change: transform;">
  <iframe loading="lazy" style="position: absolute; width: 100%; height: 100%; top: 0; left: 0; border: none; padding: 0;margin: 0;"
    src="https:&#x2F;&#x2F;www.canva.com&#x2F;design&#x2F;DAFJDj8_7PE&#x2F;view?embed" allowfullscreen="allowfullscreen" allow="fullscreen">
  </iframe>
</div>
<a href="https:&#x2F;&#x2F;www.canva.com&#x2F;design&#x2F;DAFJDj8_7PE&#x2F;view?utm_content=DAFJDj8_7PE&amp;utm_campaign=designshare&amp;utm_medium=embeds&amp;utm_source=link" target="_blank" rel="noopener">CURIS Poster</a> by Allison Tee

All major PNG-related files are in the folder labeled izhikevich_pngs. One can use the theDataMachine Python file to construct a network of neurons using the ![Brian 2](https://brian2.readthedocs.io/en/stable/) simulator and export the spike data into a suitable format for PP-(Almost-)Seq. The png and autoPng files use the PP-Almost-Seq algorithm to cluster spike data into sequences. autoPng allows the user to analyze multiple spike data files in one runthrough, and both save the resulting graph as a png (the image file format) in the graphs folder. For information on how to use autoPng and png, check out the songbird demo folder provided by the original authors.

The prototype PP-Almost-Seq code can be found in the gibbs.jl file. The author plans to clean up and sort the currently extremely long dictionary output. In the theDataMachine file, there are multiple functions to trigger sequences in data with various instances and delays. There is also a section of code that allows one to embed PNGs in the network via network connectivity and weights.

![image](https://user-images.githubusercontent.com/71671264/180329747-99292855-f9f2-4615-a699-5ce7eaad1616.png)


From the original readme:

# PP-Seq

This repo implements the point process model of neural sequences (PP-Seq) described in:

> **[Alex H. Williams](https://alexhwilliams.info) :coffee:, [Anthony Degleris](https://degleris1.github.io/) :sunrise_over_mountains:, [Yixin Wang](http://people.eecs.berkeley.edu/~ywang/), [Scott W. Linderman](https://web.stanford.edu/~swl1/) :loudspeaker:.**
> <br>[Point process models for sequence detection in high-dimensional neural spike trains](https://arxiv.org/abs/2010.04875).
> <br> *Neural Information Processing Systems 2020*, Vancouver, CA.

This model aims to identify sequential firing patterns in neural spike trains in an unsupervised manner.
For example, consider the spike train below<sup>(1)</sup>:

![image](https://user-images.githubusercontent.com/636625/94763493-36bc0400-035f-11eb-95b7-40f583ed599b.png)

By eye, we see no obvious structure in these data. However, by re-ordering the neurons according to PP-Seq's inferred sequences, we obtain:

![image](https://user-images.githubusercontent.com/636625/94767663-ee521580-0361-11eb-8729-de9eee1ab7d9.png)

Further, the model provides (probabilistic) assignment labels to each spike. In this case, we fit a model with two types of sequences. Below we use the model to color each spike as red (sequence 1), blue (sequence 2), or black (non-sequence background spike):

![image](https://user-images.githubusercontent.com/636625/94767637-db3f4580-0361-11eb-89dd-28fd8a25c468.png)

The model is fully probabilistic and Bayesian, so there are many other nice summary statistics and plots that we can make. See our paper for full details.

**Footnote.** (1) These data are deconvolved spikes from a calcium imaging recording from zebra finch HVC. These data were published in [Mackevicius*, Bahle*, et al. (2019)](https://elifesciences.org/articles/38471) and are freely available online at https://github.com/FeeLab/seqNMF.
