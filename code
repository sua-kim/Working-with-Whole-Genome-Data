%% 2020-1 바이오데이터전산기초및실습 프로젝트
% 1870021_휴먼기계바이오공학부_김수아

fprintf("============================================================\n");
fprintf("          바이오데이터전산기초및실습 기말대체프로젝트\n");
fprintf("              1870021 휴먼기계바이오공학부 김수아\n");
fprintf("============================================================\n");

%% Data Import

fprintf("============================================================\n")
fprintf("                      1. Data Import\n");
fprintf("============================================================\n")
filename = 'Homo_sapiens.GRCh38.dna.chromosome.1.fa';
% choromosome1의 정보 확인
fileInfo = dir(which(filename))

fidIn = fopen(filename, 'r');  % 파일의 포인터
header = fgetl(fidIn);
% 메모리 매핑할 파일 오픈
[fullPath, filename2, extension] = fileparts(filename);
% 매핑 파일 이름 설정
mmFilename = [filename2, '.mm'];
fidOut = fopen(mmFilename, 'w');

%% Processing & visualization

% Processing
fprintf("============================================================\n")
fprintf("                2. Processing - Memory Mapping\n");
fprintf("============================================================\n")
% 데이터를 문자 벡터로 포맷
newLine = sprintf('\n');
% 블록 사이즈 설정
blockSize = 2^20;

% end of file이 아닐 때까지 반복
while ~feof(fidIn)
    % 데이터를 문자형으로 읽어옴
    charData = fread(fidIn,blockSize,'*char')';
    % new-line 기호('\n') 제거
    charData = strrep(charData,newLine,'');
    % nt2int: 정수형으로 변환
    intData = nt2int(charData);
    % unit8로 변환후 새로운 파일에 작성
    fwrite(fidOut,intData,'uint8');
end

% 파일 닫기
fclose(fidIn);
fclose(fidOut);

% 파일을 메모리에 매핑하는 memmapfile 객체 구성 
chr1 = memmapfile(mmFilename, 'format', 'uint8');

% visualization
%
% 새롭게 생성한 파일의 정보 확인(이전과 동일하지만 newline이 제외)
mmfileInfo = dir(mmFilename)
% 인덱싱 작업을 통한 데이터 엑세스
chr1.Data(1:10)
chr1.Data(10000000:10000010)'
% 정수형으로 변환된 서열 정보를 다시 문자로 바꿔 시퀀스 정보 표시
int2nt(chr1.Data(10000000:10000010)')
seqdisp(chr1.Data(10000000:10001000)')  % 시퀀스 표시
%% Analysis 1. 블록 별 GC Content 계산하고 결과 plot
% 각 블록별 뉴클레오타이드(A,C,G,T)에 대한 GC 비율 계산

fprintf("============================================================\n")
fprintf("                         3. Analysis\n");
fprintf("============================================================\n")
fprintf("(1) 블록 별 GC Content 계산하고 결과 plot\n");
numNT = numel(chr1.Data); % 전체 길이
blockSize = 500000; % blocksize 설정
numBlocks = floor(numNT/blockSize);

% 배열을 미리 할당
ratio = zeros(numBlocks+1, 1); % 실수부분에 해당하는 데이터 처리를 위해 +1

% C와 G를 찾고 이 숫자를 A,T,C,G의 총 수로 나눔
A = nt2int('A');
C = nt2int('C');
G = nt2int('G');
T = nt2int('T');
for count = 1:numBlocks
    % block의 인덱스 계산(시작점, 끝점)
    start = 1 + blockSize*(count-1); 
    stop = blockSize*count;
    block = chr1.Data(start:stop); % block 추출
    % GC, AT 개수 찾기
    gc = (sum(block == G | block == C));
    at = (sum(block == A | block == T));
    ratio(count) = gc/(gc+at);
end

% 나머지 부분에 대한 GC-content
block = chr1.Data(stop+1: end);
gc = (sum(block == G | block == C));
at = (sum(block == A | block == T));
ratio(end) = gc/(gc+at);

% visualization: graph로 표현하기
% vector로 연결
f1 = figure;
figure(f1);
xAxis = [1:blockSize:numBlocks*blockSize, numNT];
plot(xAxis,ratio);  % 그래프 생성
xlabel('Base pairs');  % x축 라벨 설정
ylabel('Relative GC content');  % y축 라벨 설정
title('Relative GC content of Homo Sapiens Chromosome 1'); % 그래프 제목 설정

%% Analysis 2. Finding Regions of High GC Content
% GC-content가 높은 지역 찾기: GC 비율이 0.5 초과인 구간

h_indices = find(ratio>0.5);
ranges = [(1+blockSize*(h_indices-1)), blockSize*h_indices];
fprintf("\n");
fprintf("(2) Finding Regions of High GC Content\n");
fprintf('Region %d: %d has GC content %f\n', [ranges, ratio(h_indices)]');

%% Analysis 3. Finding Regions of Low GC Content
% GC-content가 낮은 지역 찾기: GC 비율이 0.4 미만인 구간

l_indices = find(ratio<0.4);
ranges = [(1+blockSize*(l_indices-1)), blockSize*l_indices];
fprintf("\n");
fprintf("(3) Finding Regions of Low GC Content\n");
fprintf('Region %d: %d has GC content %f\n', [ranges, ratio(l_indices)]');

%% Analysis 4. A,T,C,G 외의 서열 찾기
% A,T,C,G 외의 서열 찾고 블록 별 해당 서열 content 계산 뒤 결과 plot

fprintf("\n");
fprintf("(4) Analysis 5.Find sequences other than A, T, C, G\n");
no_atgc = [];  % A,T,C,G 외의 서열
for count = 1:numBlocks
    % block의 인덱스 계산(시작점, 끝점)
    start = 1 + blockSize*(count-1); 
    stop = blockSize*count;
    block = chr1.Data(start:stop); % block 추출
    % A,T,C,G 외의 서열인 경우 no_atgc 리스트에 저장, content count
    if block ~= G & block ~= C & block ~= A & block ~= T
        no_atgc = [no_atgc, block];
    end
    no_atgc_num = (sum(block ~= G & block ~= C & block ~= A & block ~= T));
    agct = (sum(block == G | block == C | block == A | block == T));
    ratio(count) = no_atgc_num/(agct+no_atgc_num);
end
% 나머지 부분에 대해서도 진행
block = chr1.Data(stop+1: end);
if block ~= G & block ~= C & block ~= A & block ~= T
    no_atgc = [no_atgc, block];
end
no_atgc_num = (sum(block ~= G & block ~= C & block ~= A & block ~= T));
agct = (sum(block == G | block == C | block == A | block == T));
ratio(end) = no_atgc_num/(agct+no_atgc_num);

%  A,T,C,G 외의 서열 출력
c = unique(no_atgc);  % 중복되는 서열 제거
fprintf("- A,T,C,G 외의 서열: %s\n", int2nt(c))

%  각 블록별 전체 서열에 대한 A,T,C,G 외의 서열 비율 그래프
fprintf("- 전체 서열에 대한 A,T,C,G외의 서열 비율 그래프")
% 앞의 Analysis 1의 그래프와 같이 출력할 수 있도록 figure 2 지정
f2 = figure;
figure(f2);
xAxis = [1:blockSize:numBlocks*blockSize, numNT];
plot(xAxis,ratio);
xlabel('Base pairs');
ylabel('Relative no AGCT content');
title('Relative no AGCT content of Homo Sapiens Chromosome 1')
