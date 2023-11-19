function train_test_ds = create_LabeledAudioDatastore(data_folder,signal_noise_ds,params)

signal_ds_main = signal_noise_ds.signal_ds_main;
noise_ds_main = signal_noise_ds.noise_ds_main;


ds_train = LabeledAudioDatastore(signal_ds_main.train,noise_ds_main.train,params);
ds_test = LabeledAudioDatastore(signal_ds_main.test,noise_ds_main.test,params);

train_test_ds.ds_train = ds_train;
train_test_ds.ds_test = ds_test;

%Save data stores

ds_file_name = sprintf('train_test_datastores_%s',datestr(now, 'yyyy_mm_dd_HH_MM'));
ds_file_name = fullfile(data_folder,ds_file_name);
save(ds_file_name,'train_test_ds');
