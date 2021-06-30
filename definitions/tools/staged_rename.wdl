version 1.0

task stagedRename {
  input {
    File original
    String name
  }

  runtime {
    memory: "4GB"
    cpu: 1
    docker: "ubuntu:bionic"
  }

  command <<<
    /bin/mv ~{original} ~{name}
  >>>

  output {
    File replacement = name
  }
}
