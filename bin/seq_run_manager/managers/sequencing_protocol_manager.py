from collections.abc import Iterator
from pathlib import Path

import minknow_api as mk
from minknow_api import protocol_pb2
from minknow_api.tools import protocols

from ..models.run_config import RunConfig


class SequencingProtocolManager:
    @staticmethod
    def get_sequencing_protocol(
        position_connection: mk.Connection,
        flowcell_product_code: str,
        kit: str,
    ) -> mk.protocol_pb2.ProtocolInfo | None:  # type: ignore
        return protocols.find_protocol(
            device_connection=position_connection,
            product_code=flowcell_product_code,
            kit=kit,
        )

    @staticmethod
    def get_available_protocols(
        position_connection: mk.Connection,
    ) -> list[mk.protocol_pb2.ProtocolInfo]:  # type: ignore
        response = position_connection.protocol.list_protocols(force_reload=True)  # type: ignore
        return response.protocols  # type: ignore

    @staticmethod
    def start_sequencing_protocol(
        device_connection: mk.Connection,
        protocol: mk.protocol_pb2.ProtocolInfo,  # type: ignore
        sample_id: str,
        sample_dir: Path,
        run_config: RunConfig,
        simplex_model: str,
        min_qscore: float,
    ) -> str:

        alignment_args = None
        if run_config.reference_genome_path is not None:
            alignment_args = protocols.AlignmentArgs(
                reference_files=[run_config.reference_genome_path],
                bed_file=run_config.sampling_regions_path,
            )
        basecalling_args = protocols.BasecallingArgs(
            simplex_model=simplex_model,
            modified_models=None,
            stereo_model=None,
            barcoding=None,
            alignment=alignment_args,
            min_qscore=min_qscore,
        )

        read_until_args = None
        if (run_config.adaptive_sampling_mode is not None) and (
            run_config.sampling_regions_path is not None
        ):
            read_until_args = protocols.ReadUntilArgs(
                filter_type=run_config.adaptive_sampling_mode,
                reference_files=[run_config.reference_genome_path],
                bed_file=run_config.sampling_regions_path,
                first_channel=None,
                last_channel=None,
            )

        pod5_args = protocols.OutputArgs(
            reads_per_file=run_config.output_chunk_size, batch_duration=None
        )
        fastq_args = protocols.OutputArgs(
            reads_per_file=run_config.output_chunk_size, batch_duration=None
        )
        bam_args = protocols.OutputArgs(reads_per_file=4000, batch_duration=None)

        protocol_identifier = protocol.identifier  # type: ignore

        protocol_args = protocols.make_protocol_arguments(
            basecalling=basecalling_args,
            read_until=read_until_args,  # type: ignore
            fastq_arguments=fastq_args,
            pod5_arguments=pod5_args,
            bam_arguments=bam_args,
            args=[],
        )

        user_info = protocol_pb2.ProtocolRunUserInfo()  # type: ignore
        user_info.protocol_group_id.value = run_config.experiment_id
        user_info.sample_id.value = sample_id

        sequencing_dir = sample_dir.parents[1]
        offload_location_info = protocol_pb2.OffloadLocationInfo()  # type: ignore
        offload_location_info.offload_location_path = sequencing_dir.as_posix()

        print(f"Offload location path: {offload_location_info.offload_location_path}")
        response = device_connection.protocol.start_protocol(  # type: ignore
            identifier=protocol_identifier,
            args=protocol_args,
            user_info=user_info,
            offload_location_info=offload_location_info,
        )
        return response.run_id

    @staticmethod
    def stream_current_protocol_updates(
        position_connection: mk.Connection,
    ) -> Iterator[mk.protocol_pb2.ProtocolRunInfo]:  # type: ignore
        return position_connection.protocol.watch_current_protocol_run()  # type: ignore

    @staticmethod
    def stream_acquisition_output(
        position_connection: mk.Connection, run_id: str
    ) -> Iterator[mk.statistics_pb2.StreamAcquisitionOutputResponse]:  # type: ignore
        return position_connection.statistics.stream_acquisition_output(acquisition_run_id=run_id)  # type: ignore

    @staticmethod
    def stop_sequencing_protocol(
        position_connection: mk.Connection, run_id: str
    ) -> protocol_pb2.StopProtocolResponse:  # type: ignore
        position_connection.protocol.stop_protocol(protocol_run_id = run_id) # type: ignore

    @staticmethod
    def get_currently_active_protocol(
        position_connection: mk.Connection,
    ) -> protocol_pb2.ProtocolRunInfo: # type: ignore
        return position_connection.protocol.get_current_protocol_run() # type: ignore
