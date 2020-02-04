% test data passing for use in web service
vector_data_in = getenv('VECTOR_DATA_IN');
vector_data_out = getenv('VECTOR_DATA_OUT');

% get input data from json file
data_in = jsondecode(fileread(vector_data_in));
disp(data_in)

% calculation using data_in array goes here
data_out = data_in;

% write output data as json
file_id = fopen(vector_data_out,'wt');
fprintf(file_id, jsonencode(data_out));
fclose(file_id);

