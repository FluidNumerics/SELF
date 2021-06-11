SELF_NAME ?= "self"
SELF_PROJECT ?= ""
SELF_ZONE ?= "us-west1-b"
SELF_MACHINE_TYPE ?= "n1-standard-2"
SELF_NODE_COUNT ?= 1
SELF_IMAGE ?= "projects/hpc-apps/global/images/openmpi-singularity-latest"
SELF_GPU_TYPE ?= "nvidia-tesla-v100"
SELF_GPU_COUNT ?= 0

.PHONY: plan apply destroy

fluid.tfvars: fluid.tfvars.tmpl
	cp fluid.tfvars.tmpl fluid.tfvars
	sed -i "s/<name>/${SELF_NAME}/g" fluid.tfvars
	sed -i "s/<project>/${SELF_PROJECT}/g" fluid.tfvars
	sed -i "s/<zone>/${SELF_ZONE}/g" fluid.tfvars
	sed -i "s/<machine_type>/${SELF_MACHINE_TYPE}/g" fluid.tfvars
	sed -i "s/<node_count>/${SELF_NODE_COUNT}/g" fluid.tfvars
	sed -i "s#<image>#${SELF_IMAGE}#g" fluid.tfvars
	sed -i "s/<gpu_type>/${SELF_GPU_TYPE}/g" fluid.tfvars
	sed -i "s/<gpu_count>/${SELF_GPU_COUNT}/g" fluid.tfvars

.terraform: 
	terraform init

plan: fluid.tfvars .terraform
	terraform plan -var-file=fluid.tfvars -out terraform.tfplan

apply: plan
	terraform apply -var-file=fluid.tfvars -auto-approve

destroy:
	terraform destroy -var-file=fluid.tfvars -auto-approve
