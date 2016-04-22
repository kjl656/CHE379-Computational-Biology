%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: BNTutorial_v2.m
% Description: Tutorial on Bayesian Networks for ChE 379 SQ16
%  This is based upon the Bayes Net Toolbox created by Kevin Murphy:
%  Reference: http://code.google.com/p/bnt/
%  Link is functional as of April 18, 2016
% Created:        April 15, 2014  by Josh Leonard
% Last modified:  April 18, 2016  by Josh Leonard
% Revision notes:
% - Update April 18, 2016: since the bioinformatics toolbox is no longer
% free for NU users, I modified the plotting functions to avoid the need
% to use the function "biograph". Hence, the plotting functions are manual
% and crude in this updated tutorial. The old code is retained below in
% case you want to use this in place of the "updated" crude code. Note that
% to use this version of the tutorial, you will need to download and load
% one additional tool, which is posted on the course website as a
% compressed file called gplotdc.zip. Unzip this file to your matlab folder
% and then add the unzipped folder to your Matlab path.
% - Update April 14, 2015: contains a fix that avoids the use of
% the BNT function "draw_graph", which is incompatible with graphics
% tools included with newer versions of MATLAB.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How to use this tutorial
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%This tutorial is divided into 4 Parts. Start by reading through each
%section before executing this script. When you run this script, MATLAB
%will pause at the end of each section. At this point, you can inspect the
%variables that were created and view their contents/values. When you are
%ready to proceed to the next section, go to the command prompt and then
%hit any key. If you would like to run this script without the pauses,
%delete the % sign before the line of code below that reads "pause off;"
%
%Note that this tutorial utilizes a random number generator. To obtain 
%repeatable results, one can set the "seed" of the random number generator
%to a fixed value, such that the results will be the same every time this
%tutorial script is run. By default, a fixed seed is *not* used, such that
%the random number generator uses a different seed each time the script is
%run. To make your results repeatable, remove the appropriate % signs as 
%indicated in Parts 2 & 3 below.

pause on;
%pause off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 1: Create a Bayesian Network
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define the topology of the graph by creating an adjacency matrix
N = 4; 
dag = false(N,N);
C=1; S=2; R=3; W=4;
dag(C,[R,S])=1;
dag(R,W)=1;
dag(S,W)=1;

%Assign names to the nodes
node_names = {'Cloudy','Sprinkler','Rain','Wetgrass'};

%Specify the size and type of each node. Here we'll let them be discrete
%and binary
discrete_nodes = 1:N;
node_sizes = [2 2 2 2]; %this says that each variable can take one of two
                        %values. This command is equivalent to writing the 
                        %following: node_sizes = 2*ones(1,N) 

%Make the Bayes Net
bnet = mk_bnet(dag, node_sizes, 'discrete', discrete_nodes);
  %this option makes bnet discrete, which is the default. We can also add
  %options to identify which nodes are "observed" if some of the variable 
  %states are unknown or "hidden" in some data sets.
  
%Now, let's draw the graph
G = bnet.dag;

%Option 1. You can draw a fancy graph using a function from the 
%bioinformatics toolbox if you have this. By default, this option is 
%commented out in the current version of this tutorial.
%gObj = biograph(G,node_names); %creates the graph object
%plot1 = view(gObj); %actually plots the graph

%Option 2. You can draw a crude manual graph using the code below. Beyond
%aestheic issues, the main disadvantage of this approach is that you must
%manually specify where to put each node on your plot using the array
%"coords" below:
coords = [2 3; 1 2; 3 2; 2 1]; %defines locations of the nodes on the plot
gplotdc(G,coords); %make the plot   
text(coords(:,1) - 0.1, coords(:,2) + 0.1, node_names(1:N)); %labels the nodes
text(0.5, 0.5, 'Nodes are connected via arcs that curve counterclockwise from start to finish');
axis([0 4 0 4]); %gives a bit more whitespace around the nodes

%Now enter model parameters by specifying the entries in the conditional
%probability table (CPD = conditional probability data).
bnet.CPD{C} = tabular_CPD(bnet, C, [0.5 0.5]);
bnet.CPD{R} = tabular_CPD(bnet, R, [0.8 0.2 0.2 0.8]);
bnet.CPD{S} = tabular_CPD(bnet, S, [0.5 0.9 0.5 0.1]);
bnet.CPD{W} = tabular_CPD(bnet, W, [1 0.1 0.1 0.01 0 0.9 0.9 0.99]);
  %If we don't specify CPD values (i.e., initialize it), the table will be
  %populated with random parameter values. If you do want to have random
  %parameters, as we'll do below, you should first initialize the random
  %number generator by first writing these 2 lines to ensure repeatability:
     %seed = 0;
     %rand('state', seed)
  %and then include lines that look like this:
     %bnet.CPD{C} = tabular_CPD(bnet, C);
     %... etc.

fprintf('\nPaused after completing Part 1. Hit any key to continue\n');
pause;
     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 2: Generate, export, and import sample data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     

%Let's now generate some sample data from this bayesian network
nsamples=100; %create 100 sample data sets from this network
samples = cell(N, nsamples);

%To obtain repeatable results, delete the % sign before the next two lines
%  seed = 0;
%  rand('state', seed);

for i=1:nsamples
  samples(:,i) = sample_bnet(bnet);
end
    %Note that we store the data as a cell array, such that each column
    %represents one "case" - a complete data set - and each row represents 
    %a node, so that each entry in the sample data table is the value that 
    %a given
    %node took in a given case.

%If you want to, you can export the sample data to a properly formatted
%text file (in ASCII format). Here, I'm first converting the samples cell
%array into a matrix of numbers (changing type of variable) so that I can
%then use dlmwrite to create my output file. Note that I've also specified
%that I want to use the space character ' ' to separate the values in my
%file since this makes the file formatted in a way that's easy to import,
%as we'll see below.
dlmwrite('sampleOutput.txt',cell2num(samples),' ');

%Now, let's load data from the text file we just created.
importedData = load('sampleOutput.txt');

fprintf('\nPaused after completing Part 2. Hit any key to continue\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 3: Parameter Estimation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Let's now see whether we can infer the parameters of a BNet using the data
%we just imported. Note that we will assume that we know the topology of
%the network, so we'll reuse the adjacency matrix that we created above
%(which was called 'dag').

%First, initialize a new bnet with random variable states:
bnet2 = mk_bnet(dag, node_sizes);

%To obtain repeatable results, delete the % sign before the next two lines
%  seed = 0;
%  rand('state', seed);

bnet2.CPD{C} = tabular_CPD(bnet2, C);
bnet2.CPD{R} = tabular_CPD(bnet2, R);
bnet2.CPD{S} = tabular_CPD(bnet2, S);
bnet2.CPD{W} = tabular_CPD(bnet2, W);

%Let's now use the maximum likelihood method to estimate the values of the
%parameters. Note that other estimation methods also exist. We'll use the
%data we just imported from our text file.
bnet3 = learn_params(bnet2, importedData);

%To view the values of the parameters that we just estimated or "learned", 
%we need to convert the CPD in bnet3 from an array of variables of type 
%tabular_CPD into a cell array of numbers. I'll do this using a little 
%'for' loop:
CPT3 = cell(1,N); %creates the cell array called CPT3
for i=1:N
  s=struct(bnet3.CPD{i});  %These commands have to do with the way that the 
                           %data objects are defined/structured in BNT
  CPT3{i}=s.CPT;
end

%Now we are ready to view the parameters. Let's look at node 4 to
%illustrate:
fprintf('\nEstimated parameters\n\n'); %prints to screen (MATLAB desktop) 
                                       %with some line breaks
dispcpt(CPT3{4});

%For comparison, let's then display the original parameters that we used
%to generate the sample data, to see how close we came with our estimate. 
%To do this, we must again convert bnet.CPD (an array of tabular_CPD 
%variables) into a cell array of numbers (i.e., a CPT): 
CPT = cell(1,N);
for i=1:N
    s=struct(bnet.CPD{i});
    CPT{i}=s.CPT;
end

%Now, we can display the original parameters to screen, preceded by a
%little caption. How close did we get?
fprintf('\nOriginal parameters\n\n');
dispcpt(CPT{4});

fprintf('\nPaused after completing Part 3. Hit any key to continue\n');
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Part 4: Topology inference
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Let's now try to infer or "learn" the topology or structure of a Bayes Net
%using our imported data set:
order = [C S R W]; %Need to specify the order of the variables in the data
                   %set that we're about to use for inference/learning
max_fan_in = 2; %Sets the largest number of parents per node. Try typing:
                %'help learn_struct_K2' at the command prompt to learn more
                % about this function.
sz = 5:5:nsamples;
for i=1:length(sz)
  dag2 = learn_struct_K2(importedData(:,1:sz(i)), node_sizes, order,... 
       max_fan_in);
end
G2 = dag2;

%Finally, let's plot the graph we just inferred/learned.

%Option 1. You can draw a fancy graph using a function from the 
%bioinformatics toolbox if you have this. By default, this option is 
%commented out in the current version of this tutorial.
%gObj2 = biograph(G2,node_names); %creates the graph object
%plot2 = view(gObj2); %actually plots the graph

%Option 2. You can draw a crude manual graph using the code below. Beyond
%aestheic issues, the main disadvantage of this approach is that you must
%manually specify where to put each node on your plot using the array
%"coords" below:

coords = [2 3; 1 2; 3 2; 2 1];
figure; %% opens a new figure window
gplotdc(G2,coords);    
text(coords(:,1) - 0.1, coords(:,2) + 0.1, node_names(1:N));
text(0.5, 0.5, 'Nodes are connected via arcs that curve counterclockwise from start to finish');
axis([0 4 0 4]);

fprintf('\nBNtutorial script complete\n');
