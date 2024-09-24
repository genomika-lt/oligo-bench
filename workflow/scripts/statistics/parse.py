import logging



logger = logging.getLogger(__name__)


def parse_sam_records(records):
    """
    Parses bam records.
    :param list[list[str]] records: Records from BAM file.
    :return: Parsed records.
    """

    parsed_records = []
    for record in records:
        meta = {}
        for meta_data in record[11:]:
            key_index = meta_data.index(':')
            type_index = meta_data.index(':', key_index + 1)
            meta_key, meta_type, meta_value = (meta_data[:key_index],
                                               meta_data[key_index + 1:type_index],
                                               meta_data[type_index + 1:])
            if meta_type == 'i':
                meta[meta_key] = int(meta_value)
            elif meta_type == 'f':
                meta[meta_key] = float(meta_value)
            elif meta_type == 'Z':
                meta[meta_key] = meta_value
            else:
                logger.warning("Unexpected Data Type: %s", meta_type)
                meta[meta_key] = meta_value

        parsed_record = record[:11] + [meta]
        parsed_records.append(parsed_record)

    return parsed_records
