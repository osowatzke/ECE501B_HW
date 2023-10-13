%% Initial Setup
% clear variables and close open figures
clear;
close all;

% set random number generator seed to get repeatable results
rng(501);

% Linear Map T
T = [1  0 -1  1  1;
     0  1  0  0  1;
    -1  0 -1 -1 -1;
     1  0 -1  0  1;
     1  1 -1  1 -1];

%% Part a)
% Create random vector in R^5
v = rand(5,1);

% Normalize vector
[~,maxIdx] = max(abs(v));
v = v/v(maxIdx);

% Apply linear map 25 times
v = applyLinearMap(T, v, 25);

%% Part b)
% Number of applications of T to consider when plotting
app = [0, 1, 3, 5, 10, 25];

% Create a figure
figure(1);
clf;
hold on;

% Plot elements of random vector for given number of applications
for i = 1:5
    plot(app,v(i,app+1),'-s','Linewidth',1.5);
end
box on;
grid on;

% Label plot
title('Value of v after Multiple Applications of T');
xlabel('Number of Applications of T')
ylabel('Value')
legendStr = cellfun(@(x) sprintf('Element %d',x), num2cell(1:5), 'UniformOutput', false);
legend(legendStr);

%% Part c)
% Create empty array for output of each run
w = zeros(5,5);

% Loop for each run
for i = 1:length(w)

    % create random vector
    v = rand(5,1);
    
    % normalize vector
    [~,maxIdx] = max(abs(v));
    v = v/v(maxIdx);
    
    % create arrays of random output v
    v = applyLinearMap(T, v, 25);

    % save last entry to array
    w(:,i) = v(:,end);
end

% Create a figure
figure(2);
clf;
hold on;

% Plot elements of random vector for each of the runs
for i = 1:5
    plot(w(i,:),'-s','Linewidth',1.5);
end
box on;
grid on;

% Label plot
title('Final Value of v for Multiple Runs');
xlabel('Run')
ylabel('Value')
legendStr = cellfun(@(x) sprintf('Element %d',x), num2cell(1:5), 'UniformOutput', false);
legend(legendStr);

%% Part d)
% Compute the eigenvalues/eigenvectors of T
[V, D] = eig(T);

% Extract the eigenvalues from the matrix D
lambda = diag(D);

% Find the index of the largest eigenvector
[~, idx] = max(abs(lambda));

% Find the eigenvector corresponding to the largest eigenvalue
w = V(:,idx);

% Normalize the resulting eigenvector
[~,maxIdx] = max(abs(w));
w = w/w(maxIdx);

%% part e)
% Reset the random number generator seed to recreate data from part a)
rng(501);

% Create random vector in R^5
v = rand(5,1);

% Normalize vector
[~,maxIdx] = max(abs(v));
v = v/v(maxIdx);

% Apply linear map 25 times
v = applyLinearMap(T, v, 25);

% Number of applications of T to consider when plotting
app = [1, 5, 10];

% Create a figure
figure(3);
clf;
hold on;

% Plot error in each element after application of T
for i = 1:5
    plot(app,abs(v(i,app+1)-w(i)),'-s','Linewidth',1.5);
end
box on;
grid on;

% Label plot
title('Error Between v and w after Repeated Applications of T');
xlabel('Number of Applications of T')
ylabel('Error')
legendStr = cellfun(@(x) sprintf('Element %d',x), num2cell(1:5), 'UniformOutput', false);
legend(legendStr);

figure(4);
clf;
hold on;

% Plot elements of random vector for given number of applications
color = cell(5,1);
for i = 1:5
    lineH = plot(app,v(i,app+1),'-s','Linewidth',1.5);
    color{i} = get(lineH, 'Color');
end

% Plot Eigenvector Corresponding to Maximum Eigenvalue for Reference
for i = 1:5
    plot([app(1) app(end)],w(i)*ones(1,2),'--','Linewidth',1.5,'Color',color{i});
end
box on;
grid on;

% Label plot
title('Comparison of v to Most Dominant Eigenvector');
xlabel('Number of Applications of T')
ylabel('Value')
legendStr = cellfun(@(x) sprintf('v(%d)',x), num2cell(1:5), 'UniformOutput', false);
legendStr2 = cellfun(@(x) sprintf('w(%d)',x), num2cell(1:5), 'UniformOutput', false);
legend([legendStr, legendStr2]);

%% Local Functions
% Function applies linear map T to vector multiple times
function v = applyLinearMap(T, v, numApplications)

    % create array of vectors input/output by T
    v = [v zeros(size(v,1),numApplications)];
    
    % Apply the linear map multiple times
    for i = 1:numApplications
    
        % Apply linear map
        v(:,i+1) = T*v(:,i);
    
        % Normalize result
        [~,maxIdx] = max(abs(v(:,i+1)));
        v(:,i+1) = v(:,i+1)/v(maxIdx,i+1);
    end
end