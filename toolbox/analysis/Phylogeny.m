%% Phylogenetic tree of Pro isolates


%% Directories
FileNames.Assembly_L1 = 'data/assemblies/Assembly_L1_20200830.mat';
FileNames.StrMod = 'data/models/StrMod.mat';
FileNames.FastaDir = 'data/genomes/IMG_Dump/nt/';

%% Load L3 Assembly, PanGEM, and StrMod
load(FileNames.Assembly_L1);
load(FileNames.StrMod);

% parse out strain names
strains = fieldnames(StrMod);
nStr = numel(strains);

%% Retrieve 16S sequence for each isolate included

strIdx = find(Pro_Assembly_L1.orgDatabase.Include);
genomeID = Pro_Assembly_L1.orgDatabase.FileID(strIdx);
nGenomes = numel(strIdx);

% number of HQ isolate genomes with 16s sequences
for a = 1:nGenomes
    strName = Pro_Assembly_L1.orgDatabase.StrainName{strIdx(a)};
    s16_idx_temp = find(strcmp('K02959',Pro_Assembly_L1.Strains.(strName).KO));
    if ~isempty(s16_idx_temp)
        s16_idx(a) = s16_idx_temp(1);
    else
        s16_idx(a) = NaN;
    end
end



for a = 1:nGenomes
    % get strain ID
    strName = Pro_Assembly_L1.orgDatabase.StrainName{strIdx(a)};
    
    % load fasta
    fastaFileName = strcat(FileNames.FastaDir,mat2str(Pro_Assembly_L1.orgDatabase.FileID(strIdx(a))),'.genes.fna');
    fasta = fastaread(fastaFileName);
    s16{a} = fasta(s16_idx(a)).Sequence;
end

