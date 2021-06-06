GMX_NAME ?= "gromacs"
GMX_PROJECT ?= ""
GMX_ZONE ?= "us-west1-b"
GMX_MACHINE_TYPE ?= "n1-standard-8"
GMX_NODE_COUNT ?= 1
GMX_IMAGE ?= "projects/hpc-apps/global/images/gromacs-gcp-dev"
GMX_GPU_TYPE ?= "nvidia-tesla-v100"
GMX_GPU_COUNT ?= 1

.PHONY: plan apply destroy

fluid.tfvars: fluid.tfvars.tmpl
	cp fluid.tfvars.tmpl fluid.tfvars
	sed -i "s/<name>/${GMX_NAME}/g" fluid.tfvars
	sed -i "s/<project>/${GMX_PROJECT}/g" fluid.tfvars
	sed -i "s/<zone>/${GMX_ZONE}/g" fluid.tfvars
	sed -i "s/<machine_type>/${GMX_MACHINE_TYPE}/g" fluid.tfvars
	sed -i "s/<node_count>/${GMX_NODE_COUNT}/g" fluid.tfvars
	sed -i "s#<image>#${GMX_IMAGE}#g" fluid.tfvars
	sed -i "s/<gpu_type>/${GMX_GPU_TYPE}/g" fluid.tfvars
	sed -i "s/<gpu_count>/${GMX_GPU_COUNT}/g" fluid.tfvars

.terraform: 
	terraform init

plan: fluid.tfvars .terraform
	terraform plan -var-file=fluid.tfvars -out terraform.tfplan

apply: plan
	terraform apply -var-file=fluid.tfvars -auto-approve

destroy:
	terraform destroy -var-file=fluid.tfvars -auto-approve
