#!/bin/bash


# Copy input deck
gsutil cp -r gs://self-fluids-data/input_decks /tmp/

cp examples/cns2d/constant_velocity/input.json /tmp/input_decks/

docker build

docker run --device=/dev/kfd \
           --device /dev/dri \
           --mount type=bind,source="/tmp/input_decks",target="/self" \
           --mount type=bind,source="$(pwd)",target="/build" \
           self:test "cd /self && /opt/self/bin/self -i /self/input.json"