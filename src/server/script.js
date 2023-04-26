
function change_command_type() {
    var use_fullcmd_checkbox = document.getElementById("use_fullcmd");
    if (use_fullcmd_checkbox.checked) {
        document.getElementById("fullcmd").disabled = false;
        document.getElementById("file").disabled = true;
        document.getElementById("slow5").disabled = true;
        document.getElementById("alignment").disabled = true;
        document.getElementById("output_dir").disabled = true;
    } else {
        document.getElementById("fullcmd").disabled = true;
        document.getElementById("file").disabled = false;
        document.getElementById("slow5").disabled = false;
        document.getElementById("alignment").disabled = false;
        document.getElementById("output_dir").disabled = false;
    }
}

function change_dir() {
    var dir_path_input = document.getElementById("dir_path");
    document.getElementById("dir_listing").src = dir_path_input.value;
}

function generate_plots() {
    document.getElementById("plot_btn").disabled = true;

    var use_fullcmd_checkbox = document.getElementById("use_fullcmd");
    var plot_command = "";
    var output_dir;
    if (use_fullcmd_checkbox.checked) {
        plot_command = document.getElementById("fullcmd").value;
        plot_cmd_parts = plot_command.split(" ")

        if (plot_cmd_parts.includes("-o")) {
            output_dir_index = plot_cmd_parts.indexOf("-o");
        } else if (plot_cmd_parts.includes("--output_dir")) {
            output_dir_index = plot_cmd_parts.indexOf("--output_dir");
        }

        output_dir = plot_cmd_parts[output_dir_index + 1];
    } else {
        var inputs = document.getElementById("plot_command_pane").getElementsByTagName("input");
        for (var i = 0; i < inputs.length; i++) {
            if (!(inputs[i].id == "use_fullcmd" || inputs[i].id == "fullcmd")) {
                plot_command += "--" + inputs[i].id + " " + inputs[i].value + " ";
            }
        }
        output_dir = document.getElementById("output_dir").value;
    }

    console.log(plot_command);

    let xhr = new XMLHttpRequest();
    xhr.open("POST", "/generate_plots");
    xhr.setRequestHeader("Accept", "application/json");
    xhr.setRequestHeader("Content-Type", "application/json");

    xhr.onreadystatechange = function () {
        if (xhr.readyState === 4) {
            document.getElementById("plot_btn").disabled = false;

            console.log(xhr.status);
            if (xhr.status == 200) {
                document.getElementById("dir_path").value = output_dir;
                document.getElementById("list_files").click();

                document.getElementById("command_history").innerHTML += "<li class=\"command-history-item\">" + plot_command + "</li>";
            } else if (xhr.status == 400) {
                alert(xhr.responseText)
            }
        }
    };

    xhr.send(plot_command);
}