# MATLAB Codes for Cycle-Free Homogeneous d-Partitions of the complete graph \(K_6\)

This repository contains the MATLAB codes used in this thesis. The codes study unordered and ordered cycle-free homogeneous \(3\)-partitions of the complete graph \(K_6\), the corresponding group actions described in Lemma 2 of Lippold et al. (2022), and elements in the groups \(H_3\) and \(G_3\).

## Mathematical Background

The computations are based on the study of cycle-free homogeneous \(d\)-partitions of \(K_{2d}\). In this work, the case \(d=3\) is considered, so the complete graph is \(K_6\).

Lemma 2 from Lippold et al. (2022) states that for a cycle-free homogeneous \(d\)-partition
\[
(\Gamma_1,\dots,\Gamma_d)
\]
of \(K_{2d}\), and for any three distinct vertices \(x,y,z\), there exists a unique cycle-free homogeneous \(d\)-partition
\[
(\Lambda_1,\dots,\Lambda_d)
\]
that agrees with the original partition on every edge except on the face \((x,y,z)\), where the two partitions differ on at least two edges.

The codes implement this operation computationally and examine the resulting permutations.

## Appendix A: Unordered Tree Partitions

### Listing 1
MATLAB code for finding unordered tree partitions of \(K_6\).

### Listing 2
MATLAB code demonstrating the group action described in Lemma 2.

### Listing 3
MATLAB code for computing the permutations in \(H_3\).

## Appendix B: Ordered Tree Partitions

### Listing 4
MATLAB code for finding ordered tree partitions of \(K_6\).

### Listing 5
MATLAB code demonstrating the group action described in Lemma 2 for ordered partitions.

### Listing 6
MATLAB code for computing the permutations in \(G_3\).

## Appendix C: Comparison of \(H_3\) and \(G_3\)

### Listing 7
MATLAB code comparing the orders of elements in \(H_3\) and \(G_3\).

## Requirements

The codes were written and tested in MATLAB. No additional MATLAB toolbox is required unless otherwise stated in the individual script files.

## How to Run

Open MATLAB and run each script individually. For example:

```matlab
run('Listing1_UnorderedPartitions_K6.m')
