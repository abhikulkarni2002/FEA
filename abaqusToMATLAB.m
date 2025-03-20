function [gcoord, nodes] = abaqusToMATLAB(inpFile)
    
    % Read the entire input file
    fileContents = fileread(inpFile);
    
    % Extract node data
    nodeData = extractBetween(fileContents, '*Node', '*Element');
    nodeLines = strtrim(splitlines(nodeData{1}));
    nodeLines = nodeLines(~cellfun('isempty', nodeLines));  % Remove empty lines
    
    % Initialize gcoord matrix
    numNodes = length(nodeLines);
    gcoord = zeros(numNodes, 2);  % Assuming 2D problem (x, y coordinates only)
    
    for i = 1:numNodes
        nodeInfo = sscanf(nodeLines{i}, '%d, %f, %f');
        if length(nodeInfo) < 3
            error('Error reading node information at line %d', i);
        end
        gcoord(i, :) = [nodeInfo(2), nodeInfo(3)];  % Store x and y coordinates
    end




    % Extract element data by searching for the *Element section
    elementStart = strfind(fileContents, '*Element');
    elementEnd = strfind(fileContents, '*Nset');
    
    if isempty(elementStart)
        error('No *Element section found in the input file.');
    end
    
    if isempty(elementEnd)
        elementEnd = length(fileContents); % If no *End Part, read until the end of the file
    end
    
    elementData = fileContents(elementStart:elementEnd); % Extract element section
    
    % Split the lines for element data
    elementLines = strtrim(splitlines(elementData));
    elementLines = elementLines(~cellfun('isempty', elementLines));  % Remove empty lines
    
    % Initialize nodes matrix for T3 elements (after splitting Q4 into 2 triangles)
    numQ4Elements = length(elementLines) - 2;  % First line is the *Element header
    nodes = zeros(numQ4Elements, 4);
    
    for i = 1:numQ4Elements
        elementInfo = sscanf(elementLines{i+1}, '%d, %d, %d, %d, %d');  % Skip the header
        
        % Make sure we have at least 5 values (ElementID, n1, n2, n3, n4)
        if length(elementInfo) < 5
            warning('Skipping invalid element line: %s', elementLines{i+1});
            continue;
        end
        
        n1 = elementInfo(2);  % First node
        n2 = elementInfo(3);  % Second node
        n3 = elementInfo(4);  % Third node
        n4 = elementInfo(5);  % Fourth node

        nodes(i, :) = [n1, n2, n3, n4];
    end