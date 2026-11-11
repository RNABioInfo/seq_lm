nextflow.enable.types = true

// The monitoring handler waits for each command before it commits state or moves on.
def sequencing_container_command(engine: String, arguments: List<String>, inputs: List<Path>,
                                 outputs: List<Path>, network: Boolean) -> List<String> {
    def image: String = 'rnabioinfo/seq_lm_seq:v1.0.0'
    def command: List<String> = []
    if (engine == 'docker') {
        def user = new com.sun.security.auth.module.UnixSystem()
        command = ['docker', 'run', '--rm', '--user', "${user.uid}:${user.gid}" as String,
                   '--group-add', '100', '--env', 'HOME=/tmp']
        command.addAll(network ? ['--add-host', 'host.docker.internal:host-gateway'] : ['--network', 'none'])
    }
    else if (engine == 'singularity' || engine == 'apptainer') {
        command = [engine, 'exec', '--cleanenv']
        image = "docker://${image}" as String
    }
    else {
        error("Sequencing monitoring does not support the '${engine}' container engine.")
    }

    inputs.toSet().each { path: Path -> command.addAll(['--mount', sequencing_bind_mount(path, true)]) }
    outputs.each { path: Path -> command.addAll(['--mount', sequencing_bind_mount(path, false)]) }
    return command + [image] + arguments
}

def sequencing_bind_mount(path: Path, read_only: Boolean) -> String {
    def source: String = path.toRealPath().toString().replace('"', '""')
    def target: String = path.toAbsolutePath().normalize().toString().replace('"', '""')
    def mode: String = read_only ? ',ro' : ''
    return "type=bind,\"src=${source}\",\"dst=${target}\"${mode}" as String
}
