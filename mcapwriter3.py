import json
import sys
import random
from time import time_ns

from mcap.writer import Writer

with open("pkpk.mcap", "wb") as stream:
    writer = Writer(stream)

    writer.start()

    schema_id = writer.register_schema(
        name="sample",
        encoding="jsonschema",
        data=json.dumps(
            {
                "type": "object",
                "properties": {
                    "sample": {
                        "type": "string",
                    }
                },
            }
        ).encode(),
    )

    channel_id = writer.register_channel(
        schema_id=schema_id,
        topic="sample_topic",
        message_encoding="json",
    )

    writer.add_message(
        channel_id=channel_id,
        log_time=time_ns(),
        data=json.dumps({"sample": '3.14'}).encode("utf-8"),
        publish_time=time_ns(),
    )
    writer.add_message(
        channel_id=channel_id,
        log_time=time_ns(),
        data=json.dumps({"sample": '4.14'}).encode("utf-8"),
        publish_time=time_ns(),
    )
    writer.add_message(
        channel_id=channel_id,
        log_time=time_ns(),
        data=json.dumps({"sample": '6.14'}).encode("utf-8"),
        publish_time=time_ns(),
    )

    writer.finish()
