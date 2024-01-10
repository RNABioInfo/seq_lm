from typing import Optional, Iterator
from minknow_api.tools import protocols
import minknow_api as mk

from models.run_config import RunConfig


class SequencingProtocolManager:
    @staticmethod
    def get_sequencing_protocol(
        position_connection: mk.Connection,
        flowcell_product_code: str,
        kit: str,
    ) -> Optional[mk.protocol_pb2.ProtocolInfo]:  # type: ignore
        return protocols.find_protocol(
            device_connection=position_connection,
            product_code=flowcell_product_code,
            kit=kit,
            basecalling=True,
        )

    @staticmethod
    def start_sequencing_protocol(
        position_connection: mk.Connection, run_config: RunConfig
    ):
        protocol = SequencingProtocolManager.get_sequencing_protocol(
            position_connection,
            run_config.flow_cell_product_code,
            run_config.kit,
        )

        if protocol is None:
            raise Exception("No protocol identifier found")

        alignment_args = protocols.AlignmentArgs(
            [run_config.reference_genome_path], run_config.sampling_regions_path
        )
        basecalling_args = protocols.BasecallingArgs(
            config=run_config.basecall_config, barcoding=None, alignment=alignment_args
        )

        read_until_args = None
        if (run_config.adaptive_sampling_mode is not None) and (
            run_config.sampling_regions_path is not None
        ):
            read_until_args = protocols.ReadUntilArgs(
                run_config.adaptive_sampling_mode,
                run_config.reference_genome_path,
                run_config.sampling_regions_path,
                None,
                None,
            )

        pod5_args = protocols.OutputArgs(run_config.output_chunk_size)
        fastq_args = protocols.OutputArgs(run_config.output_chunk_size)
        bam_args = protocols.OutputArgs(4000)

        protocol_identifier = protocol.identifier  # type: ignore

        run_id = protocols.start_protocol(
            device_connection=position_connection,
            identifier=protocol_identifier,  # type: ignore
            sample_id=run_config.sample_id,
            experiment_group=run_config.experiment_id,
            barcode_info=None,
            basecalling=basecalling_args,
            read_until=read_until_args,
            pod5_arguments=pod5_args,
            fastq_arguments=fastq_args,
            bam_arguments=bam_args,
        )

        print(f"Started run with id {run_id} on position {run_config.position_id}")

        return run_id

    @staticmethod
    def stream_current_protocol_updates(
        position_connection: mk.Connection, run_config: RunConfig
    ) -> Iterator[mk.protocol_pb2.ProtocolRunInfo]:  # type: ignore
        return position_connection.protocol.watch_current_protocol_run()  # type: ignore

    @staticmethod
    def stream_acquisition_output(
        position_connection: mk.Connection, run_id: str
    ) -> Iterator[mk.statistics_pb2.StreamAcquisitionOutputResponse]:  # type: ignore
        return position_connection.statistics.stream_acquisition_output(acquisition_run_id=run_id)  # type: ignore
