import os 
from dotenv import load_dotenv
from minio import Minio
aws_storage_bucket_name = "bio-upload-files"

path_to_input_file = os.path.join("BioinformaticsLab", "data", "data_inputs", "GenBank", 'Bacillus clausii.gb')

print(path_to_input_file)
load_dotenv()
# Initialize MinIO client
print(os.getenv("MINIO_ACCESS_ID"))
print(os.getenv("MINIO_ACCESS_PASS"))
minio_client = Minio(
    "localhost:9000",  # Replace with your MinIO server address
    access_key=os.getenv("MINIO_ACCESS_ID"),
    secret_key=os.getenv("MINIO_ACCESS_PASS"),
    secure=False  # Set to True if using HTTPS
)
# Download the object and save it to the specified file path
try:
    minio_client.fget_object(aws_storage_bucket_name, path_to_input_file, "output.gb")
    # print(f"Successfully downloaded {object_name} to {file_path}")
except Exception as e:
    print(f"Error occurred: {e}")
# client = mn.Minio(aws_storage_bucket_name, access_key=os.getenv("MINIO_ACCESS_ID"), secret_key=os.getenv("MINIO_ACCESS_PASS"), secure=True)
# minio_client.fget_object(aws_storage_bucket_name, 'output.gb', path_to_input_file)
