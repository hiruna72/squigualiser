
function change_command_type() {
    var use_fullcmd_checkbox = document.getElementById("use_fullcmd");
    var inputs = document.getElementById("plot_command_pane").getElementsByTagName("input");

    if (use_fullcmd_checkbox.checked) {
        document.getElementById("fullcmd").disabled = false;
        for (var i = 0; i < inputs.length; i++) {
            if (!(inputs[i].id == "use_fullcmd" || inputs[i].id == "fullcmd" || inputs[i].id == "use_pileup")) {
                inputs[i].disabled = true;
            }
        }
    } else {
        document.getElementById("fullcmd").disabled = true;
        for (var i = 0; i < inputs.length; i++) {
            if (!(inputs[i].id == "use_fullcmd" || inputs[i].id == "fullcmd" || inputs[i].id == "use_pileup")) {
                inputs[i].disabled = false;
            }
        }
    }
}

function change_dir() {
    var dir_path_input = document.getElementById("dir_path");
    document.getElementById("dir_listing").src = dir_path_input.value;
}

function change_dir_path(new_dir_path) {
    document.getElementById("dir_path").value = new_dir_path;
}

function generate_plots() {
    let plot_button = document.getElementById("plot_btn");
    plot_button.disabled = true;
    plot_button.innerHTML = '<span class="spinner-border spinner-border-sm" role="status" aria-hidden="true"></span> Generating...';

    var use_fullcmd_checkbox = document.getElementById("use_fullcmd");
    var use_pileup_checkbox = document.getElementById("use_pileup");
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
            if (!(inputs[i].id == "use_fullcmd" || inputs[i].id == "fullcmd" || inputs[i].id == "use_pileup")) {
                if (inputs[i].value != "") {
                    plot_command += "--" + inputs[i].id + " " + inputs[i].value + " ";
                }
            }
        }
        output_dir = document.getElementById("output_dir").value;
    }

    console.log(use_pileup_checkbox.checked + " " + plot_command);

    let xhr = new XMLHttpRequest();
    xhr.open("POST", "/generate_plots");
    xhr.setRequestHeader("Accept", "application/json");
    xhr.setRequestHeader("Content-Type", "application/json");

    xhr.onreadystatechange = function () {
        if (xhr.readyState === 4) {
            plot_button.disabled = false;
            plot_button.innerHTML = "Generate Plots";

            console.log(xhr.status);
            if (xhr.status == 200) {
                document.getElementById("error_msg").style.display = "none";
                document.getElementById("dir_listing").src = output_dir;

                document.getElementById("command_history").innerHTML = "<li>" + plot_command + "</li>" + document.getElementById("command_history").innerHTML;
            } else if (xhr.status == 400) {
                document.getElementById("error_msg").innerHTML = xhr.responseText;
                document.getElementById("error_msg").style.display = "block";
            }
        }
    };

    xhr.send(use_pileup_checkbox.checked + " " + plot_command);
}