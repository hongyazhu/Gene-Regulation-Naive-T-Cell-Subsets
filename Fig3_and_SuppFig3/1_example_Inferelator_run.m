# Example Inferelator run
# Inferelator: Miraldi et al. (2018) "Leveraging chromatin accessibility for transcriptional regulatory network inference in T Helper 17 Cells"

function funcGridSeed_lambdaTfaopt_modelSize10_priorChipKONetwork_seed()
    tfaOpts = {'', '_TFmRNA'};
    for tfaOptNum = 1:length(tfaOpts)
        for lambdaBias = [0.05, 0.25, 0.5, 0.75, 1]
            for seed = [7,99,26,57]
                rng(seed) # or rng('default') for the fifth random seed
                tfaOpt = tfaOpts{tfaOptNum};
                meanEdgesPerGene = 10;
                name_this_run = ['func_seed' num2str(seed) '_bias' num2str(100*lambdaBias) tfaOpt '_modelSize' num2str(meanEdgesPerGene)];
                diaryFile = ['outputs/230525/func_modelSize' num2str(meanEdgesPerGene) '_seed' num2str(seed) '/' name_this_run '.log'];
                matlabDir = 'infTRN_lassoStARS-master'; % downloaded from https://github.com/emiraldi/infTRN_lassoStARS
                outputFolder = ['outputs/230525/func_modelSize' num2str(meanEdgesPerGene) '_seed' num2str(seed) '/' name_this_run ];
                mkdir(['outputs/230525/func_modelSize' num2str(meanEdgesPerGene) '_seed' num2str(seed)])
                normGeneExprFile = 'gene_exp_matrix/counts_vst_blindT_combat.txt';
                targGeneFile = 'inputs/genes_DE/deg_padj01.txt';
                potRegFile = 'TF/regulators.txt';
                tfaGeneFile = 'gene_exp_matrix/counts_vst_blindT_combat.txt';
                priorName = 'chip_ko_network_atac';
                priorFile = 'prior/priors_more_than_atac_nofli1/chip_ko_network_atac.tsv';
                priorMergedTfsFile = '';
                workflow_withArguments(diaryFile, matlabDir, outputFolder, ...
                    normGeneExprFile, targGeneFile, potRegFile, tfaGeneFile, ...
                    lambdaBias, tfaOpt, priorName, priorFile, priorMergedTfsFile, ...
                    meanEdgesPerGene)
            end
        end
    end
end

%% arguments of Inferelator
% matlabDir: directory of auxiliary matlab functions (infLassoStARS, glmnet, customMatlabFxns) [..]
% outputFolder: output directory. [./outputs]
% normGeneExprFile: normalized gene expression profiles 
% targGeneFile: target genes names 
% potRegFile: potential regulators 
% tfaGeneFile: genes for TFA 
% tfaOpt: TFA estimated by equation or mRNA [ '_TFmRNA' or '' ]
% lambdaBias: prior inforcement strength [lambdaBiases = [1 .5 .25]; % correspond to no, moderate, and strong prior reinforcement]
% priorName: ['ATAC_Th17']
% priorFile: ['./inputs/priors/' priorName '.tsv']; 
% priorMergedTfsFile: ['./inputs/priors/' priorName '_mergedTfs.txt']
% meanEdgesPerGene: 10



function [] = workflow_withArguments(diaryFile, matlabDir, outputFolder, normGeneExprFile, targGeneFile, potRegFile, tfaGeneFile, lambdaBias, tfaOpt, priorName, priorFile, priorMergedTfsFile, meanEdgesPerGene, gsFile, prNickName, prTargGeneFile)

diary(diaryFile)
diary on

maxNumCompThreads(12) 

%% The following code were adapted from example_workflow_Th17 by Dr. Emily Miraldi (https://github.com/emiraldi/infTRN_lassoStARS)

%% example_workflow_Th17
% Use mLASSO-StARS to build a TRN from gene expression and prior
% information in four steps. Please refer to each function's help
% annotations for descriptions of inputs, outputs and other information.
%% References: 
% (1) Miraldi et al. (2018) "Leveraging chromatin accessibility for 
% transcriptional regulatory network inference in T Helper 17 Cells"
% (2) Qian et al. (2013) "Glmnet for Matlab."
% http://www.stanford.edu/~hastie/glmnet_matlab/
% (3) Liu, Roeder, Wasserman (2010) "Stability Approach to Regularization 
%   Selection (StARS) for High Dimensional Graphical Models". Adv. Neural.
%   Inf. Proc.
% (4) Muller, Kurtz, Bonneau. "Generalized Stability Approach for Regularized
%   Graphical Models". 23 May 2016. arXiv.
%% Author: Emily R. Miraldi, Ph.D., Divisions of Immunobiology and Biomedical
%   Informatics, Cincinnati Children's Hospital
%% Date: March 29, 2018

%clear all %%% commented the two commands!
close all
%restoredefaultpath

%argument matlabDir = '..';

addpath(fullfile(matlabDir,'infLassoStARS'))
addpath(fullfile(matlabDir,'glmnet'))
addpath(fullfile(matlabDir,'customMatlabFxns'))

%% 1. Import gene expression data, list of regulators, list of target genes
% into a Matlab .mat object
%argument geneExprTFAdir = './outputs/processedGeneExpTFA';
geneExprTFAdir = append(outputFolder, '/processedGeneExpTFA'); % added!!!

mkdir(geneExprTFAdir)
%argument normGeneExprFile = './inputs/geneExpression/th17_RNAseq254_DESeq2_VSDcounts.txt';
%argument targGeneFile = './inputs/targRegLists/targetGenes_names.txt';
%argument potRegFile = './inputs/targRegLists/potRegs_names.txt';
%tfaGeneFile = 'inputs/gene_exp_matrix/counts_vst_blindT_combat.txt';
geneExprMat = fullfile(geneExprTFAdir,'geneExprGeneLists.mat');

disp('1. importGeneExpGeneLists.m')
importGeneExpGeneLists(normGeneExprFile,targGeneFile,potRegFile,...
    tfaGeneFile,geneExprMat)

%% 2. Given a prior of TF-gene interactions, estimate transcription factor 
% activities (TFAs) using prior-based TFA and TF mRNA levels
%argument priorName = 'ATAC_Th17';
%argument priorFile = ['./inputs/priors/' priorName '.tsv']; % Th17 ATAC-seq prior
edgeSS = 50;
minTargets = 3;
[xx, priorName, ext] = fileparts(priorFile);
tfaMat = fullfile(geneExprTFAdir,[priorName '_ss' num2str(edgeSS) '.mat']);

disp('2. integratePrior_estTFA.m')
integratePrior_estTFA(geneExprMat,priorFile,edgeSS,...
     minTargets, tfaMat)

%% 3. Calculate network instabilities using bStARS

%argument lambdaBias = .5;
%argument tfaOpt = ''; % options are '_TFmRNA' or ''
totSS = 50;
targetInstability = .05;
lambdaMin = .01;
lambdaMax = 1;
extensionLimit = 1;
totLogLambdaSteps = 25; % will have this many steps per log10 within bStARS lambda range
bStarsTotSS = 5;
subsampleFrac = .63;
leaveOutSampleList = '';
leaveOutInf = '';
instabilitiesDir = fullfile(outputFolder,strrep(['instabilities_targ' ... %%% changed to outputFolder here!!!
    num2str(targetInstability) '_SS' num2str(totSS) leaveOutInf '_bS' num2str(bStarsTotSS)],'.','p'));
mkdir(instabilitiesDir)
netSummary = [priorName '_bias' strrep(num2str(100*lambdaBias),'.','p') tfaOpt];
instabOutMat = fullfile(instabilitiesDir,netSummary);

disp('3. estimateInstabilitiesTRNbStARS.m')
estimateInstabilitiesTRNbStARS(geneExprMat,tfaMat,lambdaBias,tfaOpt,...
    totSS,targetInstability,lambdaMin,lambdaMax,totLogLambdaSteps,...
    subsampleFrac,instabOutMat,leaveOutSampleList,bStarsTotSS,extensionLimit)

%% 4. For a given instability cutoff and model size, rank TF-gene
% interactions, calculate stabilities and network file for jp_gene_viz
% visualizations
%argument priorMergedTfsFile = ['./inputs/priors/' priorName '_mergedTfs.txt'];
try % not all priors have merged TFs and merged TF files
    ls(priorMergedTfsFile) 
catch
    priorMergedTfsFile = '';
end
%argument meanEdgesPerGene = 15;
targetInstability = .05;
networkDir = strrep(instabilitiesDir,'instabilities','networks');
instabSource = 'Network';
mkdir(networkDir);
networkSubDir = fullfile(networkDir,[instabSource ...
    strrep(num2str(targetInstability),'.','p') '_' ...
    num2str(meanEdgesPerGene) 'tfsPerGene']);
mkdir(networkSubDir)
trnOutMat = fullfile(networkSubDir,netSummary);
outNetFileSparse = fullfile(networkSubDir,[netSummary '_sp.tsv']);
networkHistDir = fullfile(networkSubDir,'Histograms');
mkdir(networkHistDir)
subsampHistPdf = fullfile(networkHistDir,[netSummary '_ssHist']);

disp('4. buildTRNs_mLassoStARS.m')
buildTRNs_mLassoStARS(instabOutMat,tfaMat,priorMergedTfsFile,...
    meanEdgesPerGene,targetInstability,instabSource,subsampHistPdf,trnOutMat,...
    outNetFileSparse)


%% 5. Calculate precision-recall relative to a variaty of G.S.
gsFiles = {'prior/chip_nofli1_sp.tsv', 'GS/ko_sp.tsv', 'GS/TN_network/network_sp.tsv', 'GS/tcf1_inter_sp.tsv', 'GS/tcf1_union_sp.tsv', ...
           'GS/all_tcf1Inter_sp.tsv', 'GS/all_tcf1Union_sp.tsv', 'GS/all_tcf1Inter_networks_sp.tsv', 'GS/all_tcf1Union_networks_sp.tsv'};
totGsFiles = length(gsFiles);

rankColTrn = 3;
prTargGeneFile = 'inputs/genes_DE/deg_padj01.txt';
gsRegsFile = '';
    
for gsInd = 1:totGsFiles
    gsFile = gsFiles{gsInd};
    prNickName = erase(gsFile, ["prior/", "GS/", "TN_network/", "_sp.tsv"]);
    prDir = fullfile(networkSubDir,['PR_' prNickName]);
    mkdir(prDir)
    prMatBase = fullfile(prDir,netSummary);
    prFigBase = fullfile(prDir,netSummary);

    disp('5. calcPRinfTRNs')
    calcPRinfTRNs(outNetFileSparse,gsFile,rankColTrn,...
        prTargGeneFile,gsRegsFile,prMatBase,prFigBase)
end


diary off

end
