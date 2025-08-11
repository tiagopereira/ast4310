# Projects

## Overview

There are a total of three projects, to be handed in at different times and with a different weight on the final grade, as listed below.

|        | Deadline          | Weight  |
| ------------- |:-------------:| :----:|
| Project 1  | 19 Sep | 15% |
| Project 2  | 24 Oct | 15% |
| Project 3  | 28 Nov | 20% |


## Repository

The project descriptions and questions are given in Jupyter notebook format, and can be found at the [course's github repository](https://github.com/tiagopereira/ast4310). 

Most of the time you should be using Jupyter lab or VSCode/VSCodium. In some cases the LaTeX parts in Markdown may render poorly, especially in VSCode. For Jupyter lab you can try waiting a bit, or re-running (Ctrl+Enter) the problematic markdown cells. In some cases, making a small change (e.g. adding a space) and re-running the cell fixes the problem.

## Format and delivery

The projects should be handed in as Jupyter notebooks (in Julia). The assignments will be handed using [Devilry](https://devilry.ifi.uio.no). The projects should be handed in groups of two. In exceptional cases, they can be handed individually or in groups of three.

!!! info
    Groups are not enabled in Devilry. This means that students in a group will need to submit the same files on Devilry. This is because peer review is individual (see below), therefore group members will not necessarily have the same final grade for a project.

In addition to the notebook file, you should upload also all necessary files to run a notebook, or additional files you created along your notebook (e.g. a data file with the results of a long calculation). *You do not need to upload the additional files that were present in the repository* (e.g. atom files, images, etc.). Loading data files or code from the internet is not allowed.

## Grading

### Final grade

The final grade is the 50/50 average of the final grade for all projects (numerical, from 0 to 100) and the grade from the oral exam (letter, from F to A). To convert project and oral exam grades to the final grade, we will use the following table, where the rows are different project grades and the columns different oral exam grades:

  |Project / Oral |  F |  E |  D  | C |  B |  A|
  |:---------:|----------|----|----|----|----|----|
  |**0-39**     |  F |  F |  F |  F |  F |  F |
  |**40-45**    |  F |  E |  E |  D |  D |  C |
  |**46-51**    |  F |  E |  D |  D |  C |  C |
  |**52-57**    |  F |  D |  D |  C |  C |  C |
  |**58-67**    |  F |  D |  D |  C |  C |  B |
  |**67-76**    |  F |  D |  C |  C |  B |  B |
  |**77-84**    |  F |  C |  C |  C |  B |  B |
  |**85-91**    |  F |  C |  C |  B |  B |  A |
  |**92-96**    |  F |  C |  C |  B |  B |  A |
  |**96-100**   |  F |  C |  B |  B |  A |  A |


### Project grade

The project weights are specified in the table at the top of this page. Each project will have a numerical grade between 0 and 100, and the final project grade is the weighed sum of all project grades, rounded to the closest integer. 

For the 100 points assigned to each project, 95 points are for the report itself (Jupyter notebook), and 5 points are given for peer review. 

### Peer review

To mimic the peer-review process in scientific publications, you will have the opportunity to review the work of your colleagues. The main goal of peer review is to give constructive feedback so that your colleagues can learn from their mistakes and improve the submission. **Peer review is optional.** It will give up to 5 points that will be added to your project grade. If you have a perfect score in the project but elect not to do peer review, you will get a grade of 95. 

Peer review will be double-blind and individual. Each element of a group will be given a different assignment to review. It is expected that the 5 points will be awarded for most peer reviews, although obviously low-quality peer reviews will receive zero points. Borderline cases may be given less than 5 points. Examples of a low-quality review are:

* Few comments on the work
* Mostly meaningless or off-the-mark comments, where it is obvious the reviewer did not properly read the work
* Contains offensive or inappropriate comments

Examples of meaningless comments are: *Nice work!*, *You did a good plot here!*, *Very good!*. These may be encouraging, but they don't provide any feedback to your colleagues. It's fine to give positive feedback, but please be specific, e.g.: *This part was easy to read, I enjoyed how you linked your findings to the lecture materials.*, *This is a nice plot because the two cases are clearly distinguishable, and the lines/marking are easy to see.*

#### Peer review for the *reviewer*

* After the projects are handed in in Devilry, we will work quickly to assign you a project to review. A Jupyter notebook will be sent to you on Devilry. The project deadlines are on Fridays; you should expect to receive your assignment on the following Monday.

* The deadline to submit peer review is one week after the original submission deadline (ie, Friday to Friday).

* Write your review in the notebook you receive, by adding additional markdown cells, enclosing your comments with an HTML box of the following kind (please do not change the colour of the `div` box):

```html
<div style="background-color:#e6e6e6; padding:10px; border-style:
solid;; border-color:#0d0d0d; border-width:1px">

Your markdown comments here.

* Can also use
* a bulleted list,

Other markdown, LaTeX: $\epsilon$, etc.

</div>
```

* Your job as a reviewer is to find out if each question is adequately answered, for example by answering the following questions:
    * Are the assumptions correct, and well spelled out?
    * Are the computations adequate to the job, and the code easy to read?
    * Are there any major flaws in the answer(s)?
    * Is the text well-written and concise?
    * Are plots well laid out?
* You don't need to answer all of the above, or go exhaustively through the calculations or the code. The point is more to look at the big picture.
* It is not the job of the reviewer to spell check or correct the language (in scientific publications, this is the job of the editor and/or language editor). If you find there are too many language errors, you can mention it as a general comment, e.g.: *This work would benefit from some proofreading, there are several typos and grammatical errors.*
* It is not the job of the reviewer to force her/his view into the writing of the answers, but provide a neutral viewpoint. 
* Do not be afraid to point out mistakes if you find them, but always try to do so in a constructive manner and criticise the work, not the person. Write your feedback as you would like to receive it.
* Avoid discussing your peer review assignments with your colleagues. 




#### Peer review for the *authors*

* About one week after the original submission deadline, you will receive feedback from the teachers plus any feedback from peer review
* Comments you receive from peer review will not affect your grade, only the assessment by the teachers.
* Do not take the comments from peer review at face value. It is possible that the reviewer is wrong.
* Do not take any negative comments personally. The reviewer is commenting on the work that you submitted, not on you as a person.
* Avoid discussing the comments from peer review outside your group, and especially do not confront other students or suggest she/he was the referee.

